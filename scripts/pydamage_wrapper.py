from ast import arg
import os
import subprocess
import shutil
import pandas as pd
import argparse
import time

def get_bam_files(bowtie_dir, species_under_genus):
    """Retrieve BAM files based on Genus."""
    bam_files = []
    for genus in species_under_genus["Genus"]:
        bam_file = os.path.join(bowtie_dir, f"{genus}.final.sorted.bam")
        if os.path.exists(bam_file):
            bam_files.append(bam_file)
    return bam_files

def merge_bam_files(bam_files, output_bam):
    """Merge multiple BAM files into one."""
    if not bam_files:
        print("No BAM files found!")
        return False
    merge_command = f"samtools merge -f {output_bam} {' '.join(bam_files)}"
    subprocess.run(merge_command, shell=True, check=True)
    return True

def sort_and_index_bam(input_bam, output_bam):
    """Sort and index the merged BAM file."""
    subprocess.run(f"samtools sort -o {output_bam} {input_bam}", shell=True, check=True)
    subprocess.run(f"samtools index {output_bam}", shell=True, check=True)

def run_pydamage(input_bam, output_dir):
    """Run pydamage analysis."""
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)  # Remove existing directory
    #pydamage_path = "$CONDA_PREFIX_1/envs/pydamage/bin/pydamage"
    command = f"pydamage -o {output_dir} analyze --process {process} --plot {input_bam}"
    subprocess.run(command, shell=True, check=True)

def process_pydamage_results(df, pydamage_csv):
    """Process pydamage results and update species_under_genus DataFrame."""
    if not os.path.exists(pydamage_csv):
        print(f"Error: {pydamage_csv} not found.")
        return df

    pydamage_results = pd.read_csv(pydamage_csv)
    accession_to_genus = {}

    # Create accession-to-genus mapping
    for _, row in df.iterrows():
        accessions = row["Accession_Numbers"].split(";")
        for acc in accessions:
            accession_to_genus[acc] = row["Genus"]

    df["pydamage"] = "modern"

    for _, row in pydamage_results.iterrows():
        accession = row["reference"]
        predicted_accuracy, q_val = row["predicted_accuracy"], row["qvalue"]
        genus = accession_to_genus.get(accession, None)

        # Skip if already marked as "ancient"
        if genus and df.loc[df["Genus"] == genus, "pydamage"].iloc[0] != "ancient":
            if predicted_accuracy > 0.4 and q_val <= 0.05:
                df.loc[df["Genus"] == genus, "pydamage"] = "ancient"

    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run pydamage analysis on merged BAM files.")
    parser.add_argument("species_tsv", type=str, help="Path to the species_under_genus TSV file.")
    parser.add_argument("--bowtie_dir", type=str, default="Bowtie2_results", help="Directory containing Bowtie2 BAM files.")
    parser.add_argument("--pydamage_output", type=str, default="pydamage_output", help="Directory for pydamage output.")
    parser.add_argument("--process", type=int, default=4, help="Number of process in parallel (default: 4).")

    args = parser.parse_args()
    species_tsv = args.species_tsv
    bowtie_dir = args.bowtie_dir
    pydamage_output = args.pydamage_output
    global process
    process = args.process

    if not os.path.exists(species_tsv):
        print(f"Error: The file '{species_tsv}' does not exist.")
        exit(1)

    if not os.path.exists("bowtie.finish"):
        print("Error: bowtie.finish not found. Ensure Bowtie2 processing is complete before running pydamage.")
        exit(1)

    start_time = time.time()

    species_under_genus = pd.read_csv(species_tsv, sep="\t")

    bam_files = get_bam_files(bowtie_dir, species_under_genus)
    merged_bam = "merged.bam"
    sorted_bam = "merged_sorted.bam"
    
    if merge_bam_files(bam_files, merged_bam):
        sort_and_index_bam(merged_bam, sorted_bam)
        run_pydamage(sorted_bam, pydamage_output)

        pydamage_results_path = os.path.join(pydamage_output, "pydamage_results.csv")
        species_under_genus = process_pydamage_results(species_under_genus, pydamage_results_path)
        species_under_genus.to_csv(species_tsv, sep="\t", index=False)

        # Create trigger file
        with open("pydamage.finish", "w") as f:
            f.write("pydamage analysis on test_env completed successfully.\n")

        print("\n=== pydamage Analysis Completed ===")
        print("Trigger file created: pydamage.finish")
    else:
        print("Error: Merging BAM files failed. Exiting.")
        exit(1)
        
    end_time = time.time()
    time_taken = end_time - start_time
    print(f"Total execution time: {time_taken:.2f} seconds\n\n")
