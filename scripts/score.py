import os
import pandas as pd
import argparse
import time

def check_trigger_files():
    """Ensure that all required trigger files are present before execution."""
    required_files = ["authentication.finish", "pydamage.finish", "pmdtools.finish"]
    missing_files = [file for file in required_files if not os.path.exists(file)]

    if missing_files:
        print(f"Error: Missing trigger files: {', '.join(missing_files)}")
        exit(1)

def calculate_score(row):
    """Calculate the final score based on multiple parameters"""
    score = 0
    if row["breadth_and_depth_coverage"] == "passed":
        score += 1
    if row["PMD_result"] == "ancient":
        score += 1
    if row["pydamage"] == "ancient":
        score += 1
    if row["read_dist_score"] == "ancient":
        score += 1
    score += row["edit_distance_score"]
    return score

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute final scores based on authentication parameters.")
    parser.add_argument("species_tsv", type=str, help="Path to the species_under_genus TSV file.")
    parser.add_argument("--output_file", type=str, default="final_scores.tsv", help="Output file for final scores.")
    parser.add_argument("--ancient_file", type=str, default="ancient_genera.tsv", help="Output file for ancient genera.")
    
    args = parser.parse_args()
    species_tsv = args.species_tsv
    output_file = args.output_file
    ancient_file = args.ancient_file

    if not os.path.exists(species_tsv):
        print(f"Error: The file '{species_tsv}' does not exist.")
        exit(1)

    check_trigger_files()

    start_time = time.time()

    species_under_genus = pd.read_csv(species_tsv, sep="\t")

    species_under_genus["score"] = species_under_genus.apply(calculate_score, axis=1)

    # Store the dataframe in a score-oriented format
    species_under_genus[["Genus", "species_abundance", "pydamage", "PMD_result", "breadth_and_depth_coverage",
                         "avg_read_length", "read_dist_score", "score"]].to_csv(output_file, sep="\t", index=False)
    species_under_genus.to_csv(species_tsv, sep="\t", index=False)
    ancient_genera = species_under_genus[species_under_genus["score"] > 3][["Genus", "species_abundance", "score"]]
    ancient_genera.to_csv(ancient_file, sep="\t", index=False)
    
    end_time = time.time()
    time_taken = end_time - start_time

    with open("scoring.finish", "w") as f:
        f.write("Scoring process completed successfully.\n")

    print("\n=== Final Scoring Completed ===")
    print(f"Final scores saved to: {output_file}")
    print(f"Ancient genera (final score > 3) saved to: {ancient_file}")
    print(f"Total execution time: {time_taken:.2f} seconds\n\n")

