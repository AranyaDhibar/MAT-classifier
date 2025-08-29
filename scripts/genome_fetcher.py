import os
import random
import re
import shutil
import subprocess
import zipfile
import argparse
import time
import ast
import pandas as pd 
from Bio import Entrez 
from Bio import SeqIO

global ATTEMPTS
ATTEMPTS = 0

def sanitize_name(name):
    """
    Sanitizes the species name to be safe for filenames and command-line inputs.
    Replaces special characters with underscores and removes unsafe symbols.
    """
    name = name.strip().strip('"').strip("'")
    name = re.sub(r'[\/\\\:\*\?\"<>\|]', '_', name)  # Replace unsafe characters
    name = re.sub(r'([\]\)])(?=[A-Za-z0-9])', r'\1 ', name)
    name = re.sub(r'([A-Za-z0-9])([\[\(])', r'\1 \2', name)
    name = re.sub(r'\[([^\]]+)\]', r'\1', name)  # [Genus] → Genus
    name = re.sub(r'\((.*?)\)', r'\1', name)  # (strain) → strain
    name = name.replace(" ", "_").replace(".", "")  # Remove dots and replace spaces
    name = re.sub(r'_+', '_', name)
    name = name.strip('_')
    
    return name

def normalize_name(name):
    """
    Sanitizes the species name to increase searchability.
    Omit commonly found brackets and preserve the maximum information.
    """
    name = name.strip().strip('"').strip("'")
    name = re.sub(r'([\]\)])(?=[A-Za-z0-9])', r'\1 ', name)
    name = re.sub(r'([A-Za-z0-9])([\[\(])', r'\1 \2', name)          
    name = re.sub(r'\[([^\]]+)\]', r'\1', name)         
    name = re.sub(r'\((.*?)\)', r'\1', name)              
    name = re.sub(r'\s+', ' ', name).strip()       # Normalize whitespace

    return name

def has_valid_sequences(filepath):
    """
    Check if the file has sequence and delete the file if it seems malformed
    """
    try:
        return any(len(record.seq) > 0 for record in SeqIO.parse(filepath, "fasta"))
    except Exception as e:
        os.remove(filepath)
        return False

def validate_ncbi_api_key(api_key):
    """
    Validate the provided NCBI API key.
    """
    cmd = f'datasets summary genome --limit 1 --api-key {api_key} taxon "Homo sapiens"'
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if result.returncode == 0:
        return True 

    stderr = result.stderr.lower()

    if "unauthorized" in stderr or "invalid" in stderr:
        print("Invalid API key: Unauthorized")
    else:
        print(f"API key validation failed: {stderr.strip()}")

    return False

def entrez_safe_call(callable_func, *args, **kwargs):
    ''' Runs entrez commands in a rate limited fashion'''
    global ATTEMPTS
    try:
        if not Entrez.api_key:
            time.sleep(random.uniform(0.35, 0.45))
        else:
            time.sleep(random.uniform(0.1, 0.2))
        handle = callable_func(*args, **kwargs)
        return handle

    except Exception as e:
        if any(word in str(e).lower() for word in ["server error", "gateway timeout"]):
            ATTEMPTS += 1
            print(f"NCBI service unavailable...this error occured {ATTEMPTS} times")
            if ATTEMPTS >= 15:
                print("Exiting due to repeated NCBI service unavailability. Retry later.")
                exit(1)
            
    return None

def safe_entrez_read(handle):
    """ Safely reads Entrez handle """
    try:
        return Entrez.read(handle)
    except Exception as e:
        print(f"Entrez.read() failed: {e}")
    return None

def run_datasets_command(cmd):
    """runs a shell command and returns its output"""
    global ATTEMPTS
    api_key = os.getenv('NCBI_API_KEY')
    
    try:
        if not api_key:
            time.sleep(random.uniform(0.2, 0.3))

        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode == 0 and result.stdout.strip():
            return result.stdout.strip()
        else:
            if any(word in result.stderr.strip().lower() for word in ["timeout", "temporarily", "service unavailable", "gateway timeout"]):
                ATTEMPTS += 1
                print(f"NCBI service unavailable...encountered this error {ATTEMPTS} of times")
                if ATTEMPTS >= 15:
                    print("Exiting due to repeated NCBI service unavailability. Retry later.")
                    exit(1)
    except subprocess.CalledProcessError as e:
        print(f"Command throw an error: {e}")

    return None

def extract_accession(tsv_output):
    """
    Extracts the accession number from datasets summary genome output (2nd line).
    """
    lines = tsv_output.strip().splitlines()
    if len(lines) >= 2:
        return lines[-1].strip()
    return None

def get_genome_accession(species_name):
    """
    Get the best reference genome accession for a taxon in a prioritized fallback manner
    """
    print(f"Fetching genome accession for species: {species_name}")
    cmds = [
        # Genome marked as reference: highest priority
        f'datasets summary genome --as-json-lines --limit 1 --reference taxon "{species_name}" '
        f'| dataformat tsv genome --fields accession',

        # Second prioity: RefSeq genomes 
        f'datasets summary genome --as-json-lines --limit 1 --exclude-atypical '
        f'--assembly-source "RefSeq" --assembly-level chromosome,complete taxon "{species_name}" '
        f'| dataformat tsv genome --fields accession',

        # Third priotiy: complete/chromosome assemblies but no atypical assembly
        f'datasets summary genome --as-json-lines --limit 1 --exclude-atypical  '
        f'--assembly-level chromosome,complete taxon "{species_name}" '
        f'| dataformat tsv genome --fields accession'
    ]

    for cmd in cmds:
        output = run_datasets_command(cmd)
        if output:
            acc = extract_accession(output)
            if acc and acc.startswith(("GCF_", "GCA_")):
                print(f"Genome Found for {species_name}: {acc}")
                return acc, species_name

    return search_assembly_dataset(species_name)

def search_assembly_dataset(species_name):
    """
    Searches the NCBI Assembly database for the best available genome assembly.
    """
    search_term = f"{species_name}[Organism]"

    handle = entrez_safe_call(Entrez.esearch, db="assembly", term= search_term, retmax=5, sort="relevance")
    record = safe_entrez_read(handle)
    handle.close()

    if not record["IdList"]:
        print(f"No genome assembly found for: {species_name}")
        return None

    # Fetch details of the first result
    handle = entrez_safe_call(Entrez.esummary, db="assembly", id=record["IdList"][0])
    summary_record = safe_entrez_read(handle)
    handle.close()

    # Extract assembly accession
    assembly_accession = summary_record['DocumentSummarySet']['DocumentSummary'][0].get('AssemblyAccession', 'N/A')
    
    if assembly_accession == 'N/A':
        print(f"No valid assembly accession found in 'assembly' and 'genome' database for {species_name}.")
        return None

    print(f"Genome Assembly Found in 'assembly' database for {species_name}: {assembly_accession}")
    return assembly_accession, species_name



def download_and_extract_fna(assembly_accession, species_name, output_dir="genomes", fna_output_dir="fna_files"):
    """
    Downloads and extracts the .fna genome file from NCBI RefSeq using NCBI Datasets API.
    """
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(fna_output_dir, exist_ok=True)

    # Sanitize the taxa name for filenames and paths
    safe_taxa_name = sanitize_name(species_name)
    final_fna_path = os.path.join(fna_output_dir, f"{safe_taxa_name}.fna")

    # Download the genome file
    zip_path = os.path.join(output_dir, f"{assembly_accession}.zip")
    extract_dir = os.path.join(output_dir, "extracted")
    os.makedirs(extract_dir, exist_ok=True)
    download_command = f"datasets download genome accession {assembly_accession} --include genome --filename {zip_path}"

    for attempt in range(3):
        result = subprocess.run(download_command, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error downloading {assembly_accession}:\n{result.stderr}")
            time.sleep(random.uniform(0.2, 0.3))
            continue  # Retry
        try:
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                zip_ref.extractall(extract_dir)
            break  # Success
        except zipfile.BadZipFile:
            print(f"Corrupted ZIP for {assembly_accession}, retrying (attempt {attempt + 1})...")
            if os.path.exists(zip_path):
                os.remove(zip_path)
    else:
        print(f"Failed to download valid ZIP after retries: {assembly_accession}")
        return None

    fna_file_path = None
    data_dir = os.path.join(extract_dir, "ncbi_dataset", "data", assembly_accession)

    if os.path.exists(data_dir):
        for file in os.listdir(data_dir):
            if file.endswith(".fna"):  
                fna_file_path = os.path.join(data_dir, file)
                break
            
    if fna_file_path:
        shutil.move(fna_file_path, final_fna_path)
        print(f"Extracted and renamed .fna file: {final_fna_path}")
        
        # Clean up: Remove extracted folders
        shutil.rmtree(extract_dir)
        
        return final_fna_path
    else:
        print("Error: No .fna file found in the extracted data.")
        return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch genomes from NCBI for species listed in a TSV file.")
    parser.add_argument("email", type=str, help="Required email id for NCBI Entrez record.")
    parser.add_argument("input_tsv", type=str, help="Path to the merged dataframe TSV file.")
    parser.add_argument("--api_key", type=str, default=None, help="Optional NCBI API key to increase rate limits.")
    args = parser.parse_args()
    input_tsv = args.input_tsv
    email = args.email
    api_key = args.api_key

    Entrez.email = email  

    if api_key and not validate_ncbi_api_key(api_key):
        print("WARNING: The provided NCBI API key seems invalid. Proceeding without it.")
        api_key = None
    
    if api_key:
        Entrez.api_key = api_key
        os.environ["NCBI_API_KEY"] = api_key
        
    if not os.path.exists(input_tsv):
        print(f"Error: The file '{input_tsv}' does not exist.")
        exit(1)

    start_time = time.time()

    species_under_genus = pd.read_csv(input_tsv, sep="\t")  

    required_columns = {"Genus", "scientific_name"}
    if not required_columns.issubset(species_under_genus.columns):
        print(f"Error: The input TSV must contain the columns {required_columns}")
        exit(1)

    successful_fetches = 0

    for index, row in species_under_genus.iterrows():
        genus = row["Genus"]
        species_list = row["scientific_name"]

        if isinstance(species_list, str):
            try:
                species_list = ast.literal_eval(species_list)
            except Exception as e:
                print(f"Error parsing species list for Genus '{genus}': {species_list}")
                continue

        for species_name in species_list:
            try:
                safe_taxa_name = sanitize_name(species_name)
                final_fna_path = os.path.join("fna_files", f"{safe_taxa_name}.fna")

                # Check if valid .fna file already exists
                if os.path.exists(final_fna_path):
                    if has_valid_sequences(final_fna_path):
                        successful_fetches += 1
                        continue

                species_name = normalize_name(species_name)
                accession_info = get_genome_accession(species_name)
                if accession_info:
                    assembly_accession, species_name = accession_info
                    fna_path = download_and_extract_fna(assembly_accession, species_name, output_dir="genomes", fna_output_dir="fna_files")
                    if fna_path:
                        successful_fetches += 1
            except Exception as e:
                print(f"Error processing species '{species_name}': {e}")
                continue
                
    total_attempted = sum(len(ast.literal_eval(row)) if isinstance(row, str) else len(row)
        for row in species_under_genus["scientific_name"]
    )

    end_time = time.time()
    time_taken = end_time - start_time

    print("\n=== Genome Fetching Completed ===")
    print(f"Total species attempted: {total_attempted}")
    print(f"Total {successful_fetches} genomes successfully fetched")
    print(f"Total execution time: {time_taken:.2f} seconds\n\n")