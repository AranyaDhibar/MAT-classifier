import os
import shutil
import argparse
import sys
import subprocess
import re
import pandas as pd
import threading
import psutil
import time
import datetime
import plotly.express as px

STEP_ORDER = ["genome", "concatenate", "bowtie", "pydamage", "pmdtools", "authentication", "scoring"]

def is_number_logic(value):
    try:
        float(value)
        return True
    except ValueError:
        return False
    
def script_path(script_name):
    return os.path.join(script_dir, script_name)

def expand_forced_steps(forced):
    expanded = set()
    for step in forced:
        if step in STEP_ORDER:
            index = STEP_ORDER.index(step)
            expanded.update(STEP_ORDER[index:])
    return expanded

def track_peak_resource_usage(interval=10):
    process = psutil.Process(os.getpid())
    peak = {"mem": 0.0}
    stop_event = threading.Event()

    def monitor():
        while not stop_event.is_set():
            try:
                children = process.children(recursive=True)
                all_procs = [process] + children

                total_mem = 0.0
                for p in all_procs:
                    try:
                        if p.is_running():
                            total_mem += p.memory_info().rss / (1024 ** 3)  # Convert to GB
                    except (psutil.NoSuchProcess, psutil.ZombieProcess, psutil.AccessDenied):
                        continue

                peak["mem"] = max(peak["mem"], total_mem)

            except Exception as e:
                print(f"[Monitor Warning] {e}")

            stop_event.wait(interval)

    thread = threading.Thread(target=monitor)
    thread.start()

    return stop_event, thread, peak

def parse_config(config_file):
    """Parses the config file to extract sample details and options."""
    samples = []
    options = {}

    with open(config_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):  # Ignore empty lines or comments
                continue
            
            # Detect key-value options
            if "=" in line:
                key, value = line.split("=", 1)
                key = key.strip().replace(" ", "_").lower()
                value = value.strip()
                
                # Convert numeric values if applicable
                if is_number_logic(value):  
                    value = float(value)
                    
                options[key] = value
            else:
                # Parse sample table (3 columns: sample_name, fastq, kraken_report)
                parts = line.split("\t")
                if len(parts) != 3:
                    raise ValueError(f"Error in config file:\n Expected format '<sample_name>\\t<fastq_file>\\t<kraken_report>' but got: {line}")
                
                sample_name, fastq_file, kraken_report = parts
                samples.append({
                    "sample_name": sample_name.strip(),
                    "fastq_file": fastq_file.strip(),
                    "kraken_report": kraken_report.strip()
                })

    return samples, options

def create_folders_and_copy(samples):
    """Creates folders and copies their files into them."""
    for sample in samples:
        sample_dir = sample["sample_name"]
        os.makedirs(sample_dir, exist_ok=True)

        if os.path.exists(sample["fastq_file"]):
            shutil.copy(sample["fastq_file"], sample_dir)
        else:
            print(f"Warning: Fastq file not found: {sample['fastq_file']}")

        if os.path.exists(sample["kraken_report"]):
            shutil.copy(sample["kraken_report"], sample_dir)
        else:
            print(f"Warning: Kraken report not found: {sample['kraken_report']}")

def run_kraken_parsers(samples):
    """Runs kraken_parser.py for each sample."""
    base_dir = os.getcwd()
    for sample in samples:
        sample_dir = sample["sample_name"]
        kraken_report = sample["kraken_report"]
        os.chdir(sample_dir) 
        subprocess.run(["python", script_path("kraken_parser.py"), kraken_report, "species_under_genus.tsv"], check=True)
        os.chdir(base_dir)

def merge_species_under_genus(samples):
    """Merges all species_under_genus.tsv into one grouped by Genus."""
    merged_data = {}

    for sample in samples:
        sample_dir = sample["sample_name"]
        species_file = os.path.join(sample_dir, "species_under_genus.tsv")

        if not os.path.exists(species_file):
            print(f"Warning: {species_file} not found. Skipping.")
            continue

        df = pd.read_csv(species_file, sep="\t")

        for _, row in df.iterrows():
            genus = row["Genus"]
            species_list = row["scientific_name"]

            if isinstance(species_list, str):
                species_list = species_list.strip("[]").replace("'", "").split(", ")

            if genus in merged_data:
                merged_data[genus].update(species_list)
            else:
                merged_data[genus] = set(species_list)

    merged_df = pd.DataFrame({
        "Genus": list(merged_data.keys()),
        "scientific_name": [sorted(list(species)) for species in merged_data.values()]
    })

    merged_tsv = "merged_species_under_genus.tsv"
    merged_df.to_csv(merged_tsv, sep="\t", index=False)
    print(f"\nMerged TSV saved as: {merged_tsv}")
    return merged_tsv

def run_genome_fetcher(merged_tsv, options):
    """Runs genome_fetcher.py on the merged species_under_genus.tsv file."""
    print("\nRunning Genome Fetcher on merged TSV...")
    if options.get("ncbi_api") is None:
        subprocess.run(["python", script_path("genome_fetcher.py"), options.get("email"), merged_tsv], check=True)
    else:
        subprocess.run(["python", script_path("genome_fetcher.py"), options.get("email"), merged_tsv, "--api_key", options.get("ncbi_api")], check=True)
        
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

def distribute_fna_files(samples):
    """Distributes relevant fasta files to sample directories."""
    base_fna_dir = "fna_files"

    for sample in samples:
        sample_dir = sample["sample_name"]
        species_file = os.path.join(sample_dir, "species_under_genus.tsv")
        target_fna_dir = os.path.join(sample_dir, "fna_files")
        os.makedirs(target_fna_dir, exist_ok=True)
        
        if not os.path.exists(species_file):
            print(f"Warning: {species_file} missing. Skipping {sample_dir}.")
            continue
        df = pd.read_csv(species_file, sep="\t")

        for _, row in df.iterrows():
            species_list = row["scientific_name"]
            if isinstance(species_list, str):
                species_list = species_list.strip("[]").replace("'", "").split(", ")

            for species in species_list:
                species_filename = f"{sanitize_name(species)}.fna"
                source_fna_path = os.path.join(base_fna_dir, species_filename)

                if os.path.exists(source_fna_path):
                    shutil.copy(source_fna_path, target_fna_dir)
                else:
                    print(f"Missing genome file: {species_filename}")

def run_commands(samples, options):
    """Runs the remaining commands for each sample."""
    base_dir = os.getcwd()

    for sample in samples:
        sample_dir = sample["sample_name"]
        kraken_report = sample["kraken_report"]
        fastq_file = sample["fastq_file"]
        os.chdir(sample_dir)
        print(f"\nProcessing sample: {sample_dir}...\n", flush=True)

        for step in forced_steps:
            finish_file = f"{step}.finish"
            if os.path.exists(finish_file):
                os.remove(finish_file)
                
        try:
            if not os.path.exists("concatenate.finish"):
                subprocess.run(["python", script_path("concatenate_fasta.py"), "species_under_genus.tsv", "fna_files"], check=True)
            if not os.path.exists("bowtie.finish"):
                if "bowtie" in forced_steps:
                    subprocess.run([
                    "python", script_path("bowtie_runner.py"), "species_under_genus.tsv", fastq_file,
                    "--batch_size", str(int(options.get("batch_size", 8))),
                    "--num_retries", str(int(options.get("bowtie_retries", 4))), "--force"], check=True)
                else:
                    subprocess.run([
                    "python", script_path("bowtie_runner.py"), "species_under_genus.tsv", fastq_file,
                    "--batch_size", str(int(options.get("batch_size", 8))),
                    "--num_retries", str(int(options.get("bowtie_retries", 4)))], check=True)
            if not os.path.exists("pydamage.finish"):
                subprocess.run(["python", script_path("pydamage_wrapper.py"), "species_under_genus.tsv",
                    "--process", str(int(options.get("process", 4)))], check=True)
            if not os.path.exists("pmdtools.finish"):
                subprocess.run([
                    "python", script_path("pmdtools_wrapper.py"), "species_under_genus.tsv",
                    "--pmd_threshold", str(options.get("pmd_threshold", 3.0)),  # Default: 3.0
                    "--min_fraction", str(options.get("min_fraction", 0.05)),  # Default: 0.05
                    "--pmd_library_type", options.get("pmd_library_type", "UDGhalf"), # Default: "UDGhalf"
                    "--batch_size", str(int(options.get("batch_size", 8))),  # Default: 8
                ], check=True)
            if not os.path.exists("authentication.finish"):
                subprocess.run(["python", script_path("authentication_runner.py"), "species_under_genus.tsv", kraken_report], check=True)
            if not os.path.exists("scoring.finish"):
                subprocess.run(["python", script_path("score.py"), "species_under_genus.tsv"], check=True)

            print(f"\nSuccessfully processed sample: {sample_dir}\n")
        
        except subprocess.CalledProcessError as e:
            print(f"\nError processing {sample_dir}: {e}\n")
        
        os.chdir(base_dir)

def summary(samples):
    all_samples = {}

    for sample in samples:
        sample_dir = sample["sample_name"]
        species_file = os.path.join(sample_dir, "species_under_genus.tsv")

        try:
            df = pd.read_csv(species_file, sep="\t", usecols=["Genus", "score"])
            df["score"] = pd.to_numeric(df["score"], errors="coerce")
            df.dropna(subset=["score"], inplace=True)
        except Exception as e:
            print(f"Error reading {species_file}: {e}")
            continue
        
        # Aggregate by Genus (in case there are duplicates in a sample)
        genus_scores = df.groupby("Genus")["score"].sum()

        all_samples[sample_dir] = genus_scores

    # Merge all into a single DataFrame
    merged_df = pd.DataFrame(all_samples).fillna(0).round().astype(int)
    merged_df.index.name = "Genus"
    merged_df = merged_df[(merged_df != 0).any(axis=1)]
    merged_df.to_csv("aggregated_scores.csv")
    # Filter: Drop rows where all values are ≤ 3
    merged_scores = merged_df[(merged_df > 3).any(axis=1)]
    merged_scores.to_csv('aggregated_filtered_scores.csv')

    return merged_scores
   
def summary_heatmap(merged_scores):
    if len(merged_scores) > 50:
        merged_scores["score"] = merged_scores.sum(axis=1)
        merged_df = merged_scores.nlargest(50, "score")
        title =  f"Scorewise top 50 ancient microbe candidates"
        merged_df.drop(columns="score", inplace=True)
    else:
        merged_df = merged_scores
        title = "Score distribution of putative ancient microbes"

    merged_df = merged_df.loc[merged_df.sum(axis=1).sort_values(ascending=False).index]

    fig = px.imshow(
        merged_df,
        labels=dict(x="Sample", y="Genus", color="Score"),
        x=merged_df.columns, y=merged_df.index, color_continuous_scale="YlOrRd",
        text_auto=True, aspect="auto"
    )

    fig.update_layout(
        title={'text': title, 'x': 0.5, 'xanchor': 'center', 'font': dict(size=18)},
        height=550 + len(merged_df) * 10,
        font=dict(size=12),
        coloraxis_colorbar=dict(thickness=17, title="Score", lenmode="pixels", len=200),
        xaxis=dict(tickangle=45 if len(merged_df.columns) > 5 else 0,  
            tickfont=dict(size=12), ticks="outside", tickcolor='black'),
        yaxis=dict(
            tickfont=dict(size=12 if len(merged_df) < 20 else 10),
            ticks="outside", tickcolor='black'
        )
    )
    fig.update_xaxes(automargin=True)
    fig.update_yaxes(automargin=True)
    fig.write_html(f"aggregated_putative_anc_microbes.html")
    fig.write_image(f"aggregated_putative_anc_microbes.png", scale=3)

def main():
    parser = argparse.ArgumentParser(description="Automate MAT-classifier pipeline.")
    parser.add_argument("config_file", type=str, help="Path to the configuration file.")
    parser.add_argument("--force", nargs="+", default=[], 
                        choices=["genome", "concatenate", "bowtie", "pydamage", "pmdtools", "authentication", "scoring"], 
                        help="Force specific steps to rerun. \
        Options: genome, concatenate, bowtie, pydamage, pmdtools, authentication, scoring")
    args = parser.parse_args()
    global forced_steps
    forced_steps = expand_forced_steps(set(args.force))
    
    # Start tracking resource usage
    stop_event, tracker_thread, peak_usage = track_peak_resource_usage(interval=10)
    start_time = time.time()
    try:
        samples, options = parse_config(args.config_file)
        for sample in samples:
            sample["fastq_file"] = os.path.abspath(sample["fastq_file"])
            sample["kraken_report"] = os.path.abspath(sample["kraken_report"])        
    except ValueError as e:
        print(f"\nConfiguration Error: {e}\n", file=sys.stderr)
        sys.exit(1) 
    
    print("Extracted options:")
    for key, value in options.items():
        print(f"{key} = {value}")
    
    global script_dir
    script_dir = os.path.abspath(options.get("script_dir"))
    if not os.path.isdir(script_dir):
        print(f"ERROR: script_dir path does not exist: {script_dir}")
        sys.exit(1)

    try:
        if "genome" in forced_steps or not os.path.exists("genome_fetch.finish"):
            create_folders_and_copy(samples)
            run_kraken_parsers(samples)
            merged_tsv = merge_species_under_genus(samples)
            try:
                run_genome_fetcher(merged_tsv, options)
            except subprocess.CalledProcessError as e:
                print(f"Error running genome fetcher: {e}", file=sys.stderr)
                sys.exit(1) 
            distribute_fna_files(samples)
            with open("genome_fetch.finish", "w") as f:
                f.write(f"Genome fetch and distribution completed successfully at {datetime.datetime.now()}.\n")
        forced_steps.discard("genome")

        try:
            run_commands(samples, options)
        except Exception as e:
            print(f"Unexpected error during sample processing: {e}", file=sys.stderr)
            sys.exit(1)

        try:
            summary_heatmap(summary(samples))
        except Exception as e:
            print("Error encountered during creating summary csv and plots: {e}")

        if options.get("agglomerate_and_plot") == "yes":
            agg_input_list = "all_species_under_genus_paths.txt"
            with open(agg_input_list, "w") as out_f:
                for sample in samples:
                    tsv_path = os.path.join(sample["sample_name"], "species_under_genus.tsv")
                    if os.path.exists(tsv_path):
                        out_f.write(f"{tsv_path}\n")
                    else:
                        print(f"Warning: Missing {tsv_path} — it will be skipped in agglomerate plotting.")

            if not options.get("rank") or not options.get("abundance_by"):
                print("ERROR: 'rank' and 'abundance_by' must be defined in config for agglomerate plotting.")
                print("Skipping plotting...\n")
            else:
                subprocess.run(["python", script_path("agglomerate_and_plot.py"), agg_input_list, options.get("rank"), options.get("abundance_by")], check=True)

    finally:
        # Stop the resource tracker and wait for it to finish
        stop_event.set()
        tracker_thread.join()
        end_time = time.time()
        total_time = end_time - start_time
        formatted_time = str(datetime.timedelta(seconds=int(total_time)))

        print(f"\n=== Peak Resource Usage ===")
        print(f"Total Execution Time: {formatted_time}")
        print(f"Peak Memory Usage: {peak_usage['mem']:.2f} GB")

if __name__ == "__main__":
    main()