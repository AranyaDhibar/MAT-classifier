import pandas as pd
import argparse

def parse_kraken_report(kraken_file, output_tsv):
    """
    Parses a Kraken2 output report and generates the 'species_under_genus.tsv' file.
    
    Args:
        kraken_file (str): Path to the Kraken2 report file.
        output_tsv (str): Path to save the resulting TSV file.
    """

    # Define column names for Kraken2 output
    columns = [
        "percentage", "fragments_clade", "fragments_direct", "minimizers_read",
        "distinct_minimizers", "rank_code", "taxonomy_id", "scientific_name"
    ]

    report = pd.read_csv(kraken_file, sep='\t', header=None, names=columns)

    # Filter to keep only relevant entries (clade fragments > 100 and distinct minimizers > 500)
    filtered_report = report[(report['fragments_clade'] > 100) & (report['distinct_minimizers'] > 500)].copy()
    filtered_report["scientific_name"] = filtered_report["scientific_name"].str.strip()

    genus_species = filtered_report[filtered_report["rank_code"].isin(["G", "S"])].copy()

    # Forward fill Genus names to associate species with their genera
    genus_species["Genus"] = genus_species["scientific_name"].where(genus_species["rank_code"] == "G").ffill()

    # Identify valid genera (fragments_clade > 200)
    valid_genera = report[(report["rank_code"] == "G") & (report["fragments_clade"] > 200)]["scientific_name"].str.strip().tolist()
    genus_species = genus_species[genus_species["Genus"].isin(valid_genera) | genus_species["scientific_name"].isin(valid_genera)]

    genus_species = genus_species[genus_species["rank_code"] == "S"]
    species_under_genus = genus_species.groupby("Genus")["scientific_name"].apply(list).reset_index()

    surviving_genera = set(species_under_genus["Genus"])
    eliminated_genera = set(valid_genera) - surviving_genera

    # Recover Eliminated Genera by Finding the Topmost Species under that genus
    rescued_genera = []

    for genus in eliminated_genera:
        # Find the first occurrence of the genus in the original report
        genus_row = report[(report["rank_code"] == "G") & (report["scientific_name"].str.strip() == genus)]
        
        if not genus_row.empty:
            genus_index = genus_row.index[0]  # Get the first occurrence index
            
            # Search downward for the first species (rank_code == "S")
            for i in range(genus_index + 1, len(report)):
                if report.iloc[i]["rank_code"] == "S":  # Found a species
                    rescued_genera.append({
                        "Genus": genus,
                        "scientific_name": [report.iloc[i]["scientific_name"]]
                    })
                    break  
                elif report.iloc[i]["rank_code"] == "G":  # Stop if another genus appears
                    break

    # Append recovered genera
    rescued_df = pd.DataFrame(rescued_genera)
    species_under_genus = pd.concat([species_under_genus, rescued_df], ignore_index=True)

    # Strip whitespace from species names
    species_under_genus['scientific_name'] = species_under_genus['scientific_name'].apply(
        lambda names: [name.strip() for name in names]
    )

    # Filter out rows where 'Genus' is "Homo" or ends with "virus"
    species_under_genus = species_under_genus[~species_under_genus['Genus'].str.contains(r'^(?:Homo|.*virus)$', regex=True)]
    
    species_under_genus['Genus'] = species_under_genus['Genus'].str.strip().str.replace(" ", "_")
    species_under_genus["species_abundance"] = species_under_genus["scientific_name"].apply(len)
    species_under_genus.to_csv(output_tsv, sep="\t", index=False)

    print(f"\n=== Parsing Completed for {kraken_file}===")
    print(f"Generated: {output_tsv}")
    print(f"Total Genera Processed: {len(species_under_genus)}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse Kraken2 output to generate species_under_genus.tsv.")
    parser.add_argument("kraken_file", type=str, help="Path to the Kraken2 report file.")
    parser.add_argument("output_tsv", type=str, help="Path to save the species_under_genus.tsv file.")

    args = parser.parse_args()
    
    # Run the parsing function
    parse_kraken_report(args.kraken_file, args.output_tsv)
