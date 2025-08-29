import os
import subprocess
import pandas as pd
import argparse
import time
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

def run_and_analyze_pmdtools(args):
	"""Run PMDTools and analyze score for a single genus."""
	genus, bam_file, pmd_result_dir, pmd_library_type, pmd_threshold, min_fraction = args
	output_file = os.path.join(pmd_result_dir, f"{genus}.scores")

	try:
		os.makedirs(pmd_result_dir, exist_ok=True)
		command = f"samtools view {bam_file} | pmdtools --printDS --{pmd_library_type} > {output_file}"
		subprocess.run(command, shell=True, check=True, stderr=subprocess.DEVNULL)

		if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
			df = pd.read_csv(output_file, sep="\t", header=None)
			if not df.empty:
				fraction_above_threshold = (df[3] > pmd_threshold).sum() / len(df)
				if fraction_above_threshold > min_fraction:
					return (genus, "ancient")
		return (genus, "modern")

	except subprocess.CalledProcessError:
		print(f"Error processing {genus}")
		return (genus, "unknown")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Run PMDTools analysis on BAM files.")
	parser.add_argument("species_tsv", type=str, help="Path to the species_under_genus TSV file.")
	parser.add_argument("--pmd_threshold", type=float, default=3.0, help="PMD threshold for ancient classification (default: 3.0).")
	parser.add_argument("--min_fraction", type=float, default=0.10, help="Minimum fraction of reads with PMD > threshold (default: 0.10).")
	parser.add_argument("--pmd_library_type", type=str, default="UDGhalf", help="PMDTools library type (default: UDGhalf).")
	parser.add_argument("--bowtie_dir", type=str, default="Bowtie2_results", help="Directory containing Bowtie2 BAM files.")
	parser.add_argument("--pmd_result_dir", type=str, default="pmd_results", help="Directory for PMDTools output.")
	parser.add_argument("--batch_size", type=int, default=4, help="Number of parallel PMDTools jobs.")

	args = parser.parse_args()

	if not os.path.exists(args.species_tsv):
		print(f"Error: The file '{args.species_tsv}' does not exist.")
		exit(1)
	if not os.path.exists("bowtie.finish"):
		print("Error: Ensure bowtie2 is successfully processed before running PMDTools.")
		exit(1)

	start_time = time.time()

	species_under_genus = pd.read_csv(args.species_tsv, sep="\t")
	species_under_genus["species_abundance"] = pd.to_numeric(
		species_under_genus["species_abundance"], errors="coerce").fillna(1).clip(lower=1)

	# Create args list
	pool_args = []
	for _, row in species_under_genus.iterrows():
		genus = row["Genus"]
		abundance = row["species_abundance"]

		bam_path = os.path.join(args.bowtie_dir, f"{genus}.final.sorted.bam")
		if os.path.exists(bam_path):
			# Dynamically adjust min_fraction based on abundance
			dynamic_min_fraction = max(0.01, args.min_fraction / abundance)
			
			pool_args.append((genus, bam_path, args.pmd_result_dir, args.pmd_library_type,
				args.pmd_threshold, dynamic_min_fraction))


	results = {}
	with ProcessPoolExecutor(max_workers=args.batch_size) as executor:
		for genus, result in tqdm(executor.map(run_and_analyze_pmdtools, pool_args), total=len(pool_args)):
			results[genus] = result

	species_under_genus["PMD_result"] = species_under_genus["Genus"].map(lambda g: results.get(g, "unknown"))
	species_under_genus.to_csv(args.species_tsv, sep="\t", index=False)

	with open("pmdtools.finish", "w") as f:
		f.write("PMDTools analysis completed successfully.\n")

	print("\n=== PMDTools Analysis Completed ===")
	print(f"Total execution time: {time.time() - start_time:.2f} seconds\n")
	print("Trigger file created: pmdtools.finish \n\n")