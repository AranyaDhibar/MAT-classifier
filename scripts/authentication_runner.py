import os
import subprocess
import pandas as pd
import numpy as np
import time
import pysam
from collections import Counter
from sklearn.cluster import KMeans
import argparse

def filter_kraken2_report(kraken_report_path):
    """Filter Kraken2 report based on predefined thresholds. [need atleast 200 reads & 1000 unique kmers]"""
    if not os.path.exists(kraken_report_path):
        print(f"Error: Kraken2 report '{kraken_report_path}' not found.")
        exit(1)

    columns = [
        "percentage", "fragments_clade", "fragments_direct", "minimizers_read",
        "distinct_minimizers", "rank_code", "taxonomy_id", "scientific_name"
    ]

    report = pd.read_csv(kraken_report_path, sep="\t", header=None, names=columns)
    filtered_report = report[(report['fragments_clade'] > 200) & (report['distinct_minimizers'] > 1000)].copy()
    filtered_genera = set(filtered_report[filtered_report['rank_code'] == 'G']['scientific_name'].str.strip())

    return filtered_genera

def get_bam_files(bowtie_dir, species_under_genus):
    """Retrieve BAM files for each genus."""
    bam_files = {}
    for genus in species_under_genus["Genus"]:
        bam_file = os.path.join(bowtie_dir, f"{genus}.final.sorted.bam")
        if os.path.exists(bam_file):
            bam_files[genus] = bam_file
    return bam_files

def get_avg_read_length(bam_file):
    """Calculate the average read length from a BAM file using samtools."""
    try:
        cmd = f"samtools view {bam_file} | awk '{{if ($10 == \"*\") {{sum+=0; count+=0}} else {{sum+=length($10); count+=1}}}} END {{if (count > 0) print sum/count; else print 0}}'"
        avg_length = subprocess.check_output(cmd, shell=True, text=True).strip()
        return float(avg_length) if avg_length else 0
    except Exception as e:
        return None

def perform_kmeans_clustering(species_under_genus):
    """Perform K-means clustering to classify ancient vs modern based on read length."""
    df = species_under_genus.copy()
    
    # Identify valid read lengths (non-zero & non-null)
    valid_mask = (df['avg_read_length'].notna()) & (df['avg_read_length'] > 0)
    valid_read_lengths = df.loc[valid_mask, 'avg_read_length'].values.reshape(-1, 1)
    if valid_read_lengths.shape[0] < 2:
        # If there's not enough valid data to cluster
        df['read_dist_score'] = 'unknown'
        return df

    # Calculate interquartile range to see if there is sufficient variation
    q1 = np.percentile(valid_read_lengths, 25)
    q3 = np.percentile(valid_read_lengths, 75)
    iqr = q3 - q1

    if iqr <= 3:    #3 is chosen arbitarily 
        df.loc[valid_mask, 'read_dist_score'] = 'modern'
        df.loc[~valid_mask, 'read_dist_score'] = 'unknown'
        return df

    # Apply K-means with 2 clusters
    kmeans = KMeans(n_clusters=2, random_state=32, n_init=50)
    clusters = kmeans.fit_predict(valid_read_lengths)
    df.loc[valid_mask, 'cluster_label'] = clusters

    # Identify the cluster with the smallest average read length
    smallest_cluster_label = np.argmin(kmeans.cluster_centers_.flatten())

    # Assign classification labels
    df['read_dist_score'] = df['cluster_label'].apply(
        lambda cluster_index: 'ancient' if cluster_index == smallest_cluster_label else 'modern'
    )
    df.loc[~valid_mask, 'read_dist_score'] = 'unknown'
    df.drop(columns=['cluster_label'], inplace=True)
    return df

def edit_distance_score(genus, bowtie_dir):
    """Compute edit distance score from BAM file."""
    try:
        bam_file = os.path.join(bowtie_dir, f"{genus}.final.sorted.bam")
        samfile = pysam.AlignmentFile(bam_file, "rb")
        edit_distance_counts = Counter()

        for read in samfile:
            if read.has_tag("NM"):
                edit_distance = read.get_tag("NM")
                if edit_distance in [0, 1, 2, 3, 4]: 
                    edit_distance_counts[edit_distance] += 1

        # Assign a score based on edit distance pattern
        if (edit_distance_counts.get(1, 0) > edit_distance_counts.get(2, 0) > 
            edit_distance_counts.get(3, 0) > edit_distance_counts.get(4, 0)):
            return 1
        else: 
            return 0  
    except Exception as e:
        return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute authentication parameters for genomic data.")
    parser.add_argument("species_tsv", type=str, help="Path to the species_under_genus TSV file.")
    parser.add_argument("kraken_report", type=str, help="Path to the Kraken2 report file.")
    parser.add_argument("--bowtie_dir", type=str, default="Bowtie2_results", help="Directory containing Bowtie2 BAM files.")

    args = parser.parse_args()
    species_tsv = args.species_tsv
    kraken_report = args.kraken_report
    bowtie_dir = args.bowtie_dir

    if not os.path.exists(species_tsv):
        print(f"Error: The file '{species_tsv}' does not exist.")
        exit(1)

    if not os.path.exists("bowtie.finish"):
        print("Error: Ensure bowtie2 is successfully processed before running authentication.")
        exit(1)

    start_time = time.time()
    species_under_genus = pd.read_csv(species_tsv, sep="\t")

    # Auth 1: Filter Kraken2 report and assign breadth_and_depth_coverage
    filtered_genera = filter_kraken2_report(kraken_report)
    species_under_genus['breadth_and_depth_coverage'] = species_under_genus['Genus'].apply(
        lambda x: 'passed' if x in filtered_genera else 'failed'
    )

    # Auth 2: Bin the read length into two cluster to score modern or ancient
    bam_files = get_bam_files(bowtie_dir, species_under_genus)
    species_under_genus['avg_read_length'] = species_under_genus['Genus'].apply(
        lambda x: get_avg_read_length(bam_files.get(x, '')) if x in bam_files else None
    )

    species_under_genus = perform_kmeans_clustering(species_under_genus)

    # Auth 3: Compute edit distance scores
    species_under_genus["edit_distance_score"] = species_under_genus["Genus"].apply(
        lambda x: edit_distance_score(x, bowtie_dir))

    species_under_genus.to_csv(species_tsv, sep="\t", index=False)

    end_time = time.time()
    time_taken = end_time - start_time

    with open("authentication.finish", "w") as f:
        f.write("Authentication parameters computation completed successfully.\n")

    print("\n=== Authentication Parameters Computation Completed ===")
    print(f"Total execution time: {time_taken:.2f} seconds")
    print("Trigger file created: authentication.finish\n\n")
