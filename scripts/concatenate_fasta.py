import os
import pandas as pd
import argparse
import time
import ast
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def read_and_concatenate(fasta_files):
    """
    Reads multiple FASTA files and concatenates all sequences into one.
    Also extracts the accession numbers from the headers.
    """
    concatenated_sequence = Seq("")
    first_header = None
    accession_numbers = []

    for fasta_file in fasta_files:
        with open(fasta_file, "r") as f:
            records = list(SeqIO.parse(f, "fasta"))
            if records:
                if first_header is None:
                    first_header = records[0].id  # Save the first header
                    accession_numbers.extend([record.id for record in records])
                concatenated_sequence += Seq("".join(str(record.seq) for record in records))

    if first_header is None:
        first_header = "concatenated_genome"

    return first_header, concatenated_sequence, accession_numbers

def concatenate_genus_fasta(genus, fasta_files, output_dir="fna_files"):
    """
    Concatenates multiple FASTA files for a genus and writes the output.
    """
    os.makedirs(output_dir, exist_ok=True)

    if len(fasta_files) == 1:
        # Single genome, use its file directly
        output_file = fasta_files[0]
        with open(output_file, "r") as f:
            records = list(SeqIO.parse(f, "fasta"))
            accession_numbers = [record.id for record in records]
        return output_file, accession_numbers, len(records[0].seq) if records else 0

    else:
        # Read and concatenate all sequences
        header, concatenated_seq, accession_numbers = read_and_concatenate(fasta_files)

        output_file = os.path.join(output_dir, f"{genus}_concatenated.fna")
        new_record = SeqRecord(concatenated_seq, id=header, description=f"Concatenated sequence for {genus}")

        # Write to output file
        with open(output_file, "w") as out_f:
            SeqIO.write(new_record, out_f, "fasta")

        return output_file, accession_numbers, len(concatenated_seq)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Concatenate FASTA files for genera with multiple species.")
    parser.add_argument("species_tsv", type=str, help="Path to the species_under_genus TSV file.")
    parser.add_argument("genome_dir", type=str, help="Path to the directory containing genome FASTA files.")

    args = parser.parse_args()
    species_tsv = args.species_tsv
    genome_dir = args.genome_dir

    if not os.path.exists(species_tsv):
        print(f"Error: The file '{species_tsv}' does not exist.")
        exit(1)

    if not os.path.exists(genome_dir):
        print(f"Error: The genome directory '{genome_dir}' does not exist.")
        exit(1)

    start_time = time.time() 

    species_under_genus = pd.read_csv(species_tsv, sep="\t")

    results = []

    for index, row in species_under_genus.iterrows():
        genus = row["Genus"]
        species_list = ast.literal_eval(row["scientific_name"]) 

        fasta_files = []
        
        # Find FASTA files for the given species names
        for species_name in species_list:
            safe_name = species_name.replace(" ", "_")  # Match filename convention
            fasta_path = os.path.join(genome_dir, f"{safe_name}.fna")
            
            if os.path.exists(fasta_path):
                fasta_files.append(fasta_path)
        
        # Concatenate files for the genus
        output_file, accession_numbers, sequence_length = concatenate_genus_fasta(genus, fasta_files, output_dir=genome_dir)

        # Store results
        results.append({
            "Genus": genus,
            "Output_File": output_file,
            "Accession_Numbers": ";".join(accession_numbers) if accession_numbers else "",
            "Sequence_Length": sequence_length
        })

    # Convert results to a DataFrame and merge with original df
    df_results = pd.DataFrame(results)
    species_under_genus = species_under_genus.merge(df_results, on="Genus", how="left")

    # Remove genera that has no valid sequence
    species_under_genus = species_under_genus[species_under_genus['Sequence_Length'] != 0]
    species_under_genus.to_csv(species_tsv, sep="\t", index=False)

    end_time = time.time()

    time_taken = end_time - start_time

    # Create a trigger file
    trigger_file = "concatenate.finish"
    with open(trigger_file, "w") as f:
        f.write("Concatenation process completed successfully.\n")

    print("\n=== FASTA Concatenation Completed ===")
    print(f"Total Genera Processed: {len(species_under_genus)}")
    print(f"Total execution time: {time_taken:.2f} seconds\n\n")

