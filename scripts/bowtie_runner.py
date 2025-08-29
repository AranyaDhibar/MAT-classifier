import os
import glob
import subprocess
import pandas as pd
import argparse
import time
from colorama import Fore, init
from concurrent.futures import ProcessPoolExecutor

def process_species(args):
    """Process a single reference genome using Bowtie2."""
    row, merged_reads, output_dir, log_dir, force = args
    if isinstance(row, dict):
        row = pd.Series(row)
    genus = row["Genus"] 
    ref_genome_prefix = row["Output_File"]  

    final_sorted_bam = f"{output_dir}/{genus}.final.sorted.bam"
    bowtie_log = os.path.join(log_dir, f"{genus}.log")
    
    if not force and os.path.exists(f"{final_sorted_bam}.bai"):
        #print(f"Skipping {genus}, already processed.")
        return genus
    
    try:
        # Step 1: Index the reference genome
        subprocess.run(
            f"bowtie2-build {ref_genome_prefix} {ref_genome_prefix}",
            shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True
        )
        
        # Step 2: Align reads using Bowtie2
        sam_file = f"{output_dir}/{genus}.sam"
        with open(bowtie_log, "w") as log_f:
            subprocess.run(
                f"bowtie2 -x {ref_genome_prefix} --end-to-end --threads 4 --very-sensitive "
                f"-U {merged_reads} -S {sam_file}",
                shell=True, stderr=log_f, check=True
            )
        
        # Step 3: Remove duplicate headers and clean SAM
        clean_sam = f"{output_dir}/{genus}.nodups.sam"
        subprocess.run(
            f"grep @ {sam_file} | awk '!seen[$2]++' > {output_dir}/{genus}_header_nodups.txt && "
            f"grep -v '^@' {sam_file} > {output_dir}/{genus}_noheader.sam && "
            f"cat {output_dir}/{genus}_header_nodups.txt {output_dir}/{genus}_noheader.sam > {clean_sam}",
            shell=True, check=True
        )

        # Step 4: Convert to BAM and sort
        subprocess.run(
            f"samtools view -bS -q 1 -h {clean_sam} | samtools sort -@ 4 -o {final_sorted_bam}",
            shell=True, check=True
        )
        subprocess.run(f"samtools index {final_sorted_bam}", shell=True, check=True)

        # Cleanup intermediate files
        os.remove(sam_file)
        os.remove(clean_sam)
        os.remove(f"{output_dir}/{genus}_header_nodups.txt")
        os.remove(f"{output_dir}/{genus}_noheader.sam")
    
    except subprocess.CalledProcessError as e:
        print(f"Error processing {genus}: {e}")
    
    return genus

def merge_logs(log_dir, merged_log_file):
    """Merge all individual log files into a single log file."""
    log_files = sorted(glob.glob(os.path.join(log_dir, "*.log")))
    with open(merged_log_file, "w") as outfile:
        for log_file in log_files:
            with open(log_file, "r") as infile:
                outfile.write(f"\n===== LOG FILE: {os.path.basename(log_file)} =====\n")
                outfile.write(infile.read())
                outfile.write("\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Bowtie2 alignment pipeline.")
    parser.add_argument("species_tsv", type=str, help="Path to the species_under_genus TSV file.")
    parser.add_argument("merged_reads", type=str, help="Path to the merged FASTQ reads file.")
    parser.add_argument("--batch_size", type=int, default=4, help="Number of species to process in parallel (default: 4).")
    parser.add_argument("--num_retries", type=int, default=3, help="Number of reruns if alignment fails (default: 3).")
    parser.add_argument("--force", action="store_true", help="Force reprocessing even if BAM index files exist.")

    args = parser.parse_args()
    species_tsv = args.species_tsv
    merged_reads = args.merged_reads
    batch_size = args.batch_size
    num_retries = args.num_retries
    num_jobs = min(batch_size, os.cpu_count() // 4)

    if not os.path.exists(species_tsv) or not os.path.exists(merged_reads):
        print(f"Error: Required input file(s) not found.")
        exit(1)

    if not os.path.exists("concatenate.finish"):
        print("Error: concatenate.finish not found. Ensure concatenation step is completed first.")
        exit(1)

    # Start time and memory tracking
    start_time = time.time()

    # Output directories
    output_dir = "Bowtie2_results"
    os.makedirs(output_dir, exist_ok=True)
    log_dir = os.path.join(output_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)
    merged_log_file = os.path.join(output_dir, "alignment.log")
    merged_bowtie_log_file = os.path.join(output_dir, "alignment_results.log")


    species_under_genus = pd.read_csv(species_tsv, sep="\t")
    expected_files = len(species_under_genus)
    
    for attempt in range(num_retries):
        species_list = [row.to_dict() for _, row in species_under_genus.iterrows()]
        args_list = [(row, merged_reads, output_dir, log_dir, args.force) for row in species_list]
        with ProcessPoolExecutor(max_workers=num_jobs) as executor:
            results = list(executor.map(process_species, args_list))

        # Batch processing with threading
        # species_list = [row for _, row in species_under_genus.iterrows()]
        # batches = [species_list[i:i + batch_size] for i in range(0, len(species_list), batch_size)]
        # with ThreadPoolExecutor(max_workers=8) as executor:
        #     results = list(executor.map(lambda batch: [process_species(row) for row in batch], batches))

        # Count completed .bai files
        bai_files = len(glob.glob(os.path.join(output_dir, "*.bai")))

        if bai_files == expected_files:
            print(f"All {bai_files}/{expected_files} alignments completed successfully.")
            break
        else:
            print(f"Retry {attempt+1}/{num_retries}: Only {bai_files}/{expected_files} alignments completed.")
    
    # If retries are exhausted, filter out failed taxa
    bai_files = len(glob.glob(os.path.join(output_dir, "*.bai")))
    if bai_files != expected_files:
        print(f"Warning: Alignment failed for some species even after {num_retries} retries.")
        completed_species = set([os.path.basename(f).split(".")[0] for f in glob.glob(os.path.join(output_dir, "*.bai"))])
        species_under_genus = species_under_genus[species_under_genus["Genus"].isin(completed_species)]
        species_under_genus.to_csv(species_tsv, sep="\t", index=False)

    # Merge logs
    merge_logs(log_dir, merged_bowtie_log_file)

    end_time = time.time()
    time_taken = end_time - start_time

    init(autoreset=True)    

    if bai_files != expected_files:
        print(Fore.RED + f"Warning: Alignment failed for some species even after {num_retries} retries.")
        completed_species = set([os.path.basename(f).split(".")[0] for f in glob.glob(os.path.join(output_dir, "*.bai"))])
        failed_species = set(species_under_genus["Genus"]) - completed_species
        print("\n=== FAILED ALIGNMENTS ===")
        for species in failed_species:
            print(f"- {species}")
        print("\nUpdated species_under_genus.tsv saved with only successfully aligned species.")

    with open("bowtie.finish", "w") as f:
        f.write("Bowtie2 alignment completed successfully!\n")

    print("\n=== Bowtie2 Alignment Completed ===")
    print(f"Total Genera Processed: {bai_files}/{expected_files}")
    print(f"Total execution time: {time_taken:.2f} seconds\n\n")
    
