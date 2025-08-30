import argparse
import ast
import os
import pandas as pd
import plotly.express as px
from Bio import Entrez
import time
from ete3 import NCBITaxa

Entrez.email = "aranyad@utexas.edu"
ncbi = NCBITaxa() 

def rank_name1(species, rank):
    taxid = ncbi.get_name_translator([species]).get(species)
    if not taxid:
        return rank_name2(species, rank)

    lineage = ncbi.get_lineage(taxid[0])
    names = ncbi.get_taxid_translator(lineage)
    ranks = ncbi.get_rank(lineage)

    for tid in lineage:
        if ranks[tid] == rank:
            return names[tid]

    return rank_name2(species, rank)

def rank_name2(species, rank):
    time.sleep(0.34)
    query = f"{species}[Scientific Name]"
    handle = Entrez.esearch(db="taxonomy", term=query)
    record = Entrez.read(handle)
    if not record["IdList"]:
        return "unknown"

    taxid = record["IdList"][0]
    handle = Entrez.efetch(db="taxonomy", id=taxid)
    result = Entrez.read(handle)
    lineage = result[0]["LineageEx"]
    for entry in lineage:
        if entry["Rank"] == rank:
            return entry["ScientificName"]

    return "unknown"

def grouping_by_abundance(df, rank, abundance_by):
    if abundance_by == "species":
        grouped_df = df.groupby(rank)["species_abundance"].sum().reset_index()
        grouped_df.rename(columns={"species_abundance": "abundance_value"}, inplace=True)
    elif abundance_by == "genus":
        grouped_df = df.groupby(rank)["Genus"].nunique().reset_index()
        grouped_df.rename(columns={"Genus": "abundance_value"}, inplace=True)

    return grouped_df

def merge_dfs(dfs, rank):
    merged_df = None
    for key, df in dfs.items():
        df.rename(columns={"abundance_value": key}, inplace=True)
        if merged_df is None:
            merged_df = df.copy()
        else:
            merged_df = pd.merge(merged_df, df, on=rank, how='outer')

    merged_df.fillna(0, inplace=True)

    return merged_df

def heatmap_plotly(df, output_name, rank, abundance_by):
    merge_df = df.copy()
    merge_df.set_index(rank, inplace=True)
    if len(merge_df) > 50:
        merge_df["total_abundance"] = merge_df.sum(axis=1)
        merge_df = merge_df.nlargest(50, "total_abundance")
        title =  f"Samplewise top 50 {abundance_by} abundance at {rank} level"
        merge_df.drop(columns="total_abundance", inplace=True)
    else:
        title =  f"Samplewise {abundance_by} abundance at {rank} level"

    fig = px.imshow(
        merge_df,
        labels=dict(x="Sample", y=rank.capitalize(), color="Abundance"),
        x=merge_df.columns, y=merge_df.index, color_continuous_scale="inferno_r",
        text_auto=True, aspect="auto"
    )

    fig.update_layout(
        title={'text': title, 'x': 0.5, 'xanchor': 'center', 'font': dict(size=18)},
        height=550 + len(merge_df) * 12,
        font=dict(size=12),
        coloraxis_colorbar=dict(thickness=17, title="Abundance", lenmode="pixels", len=200),
        xaxis=dict(tickangle=45 if len(merge_df.columns) > 5 else 0,  
            tickfont=dict(size=12), ticks="outside", tickcolor='black'),
        yaxis=dict(
            tickfont=dict(size=12 if len(merge_df) < 20 else 10),
            ticks="outside", tickcolor='black'
        )
    )

    fig.update_xaxes(automargin=True)
    fig.update_yaxes(automargin=True)
    fig.write_html(f"{output_name}.html")
    fig.write_image(f"{output_name}.png", scale=3)

def barplot_plotly(df, output_name, rank, abundance_by):
    # taking 10 top ranks and plot other ranks as "others"
    df = df.set_index(rank)
    abundances = df.sum(axis=1).sort_values(ascending=False)
    top_taxa = abundances.head(10).index 
    df_trimmed = df.copy()
    df_trimmed['Taxon'] = df_trimmed.index
    df_trimmed['Taxon'] = df_trimmed['Taxon'].apply(lambda x: x if x in top_taxa else 'Others')
    df_long = df_trimmed.melt(id_vars='Taxon', var_name='Sample', value_name='Abundance')

    fig = px.bar(
        df_long, x='Sample', y='Abundance', color='Taxon',
        color_discrete_sequence=px.colors.qualitative.Pastel,
        category_orders={'Taxon': list(top_taxa) + ['Others']},
        title=f'{rank.capitalize()} Abundance per Sample',
        labels={'Abundance': f'{abundance_by.capitalize()} Abundance', 'Sample': 'Sample', 'Taxon': rank.capitalize()},
    )

    fig.update_layout(
        title={'x': 0.5, 'xanchor': 'center', 'font': dict(size=18)},
        legend=dict(
            orientation="h", yanchor="bottom", y=-0.4, xanchor="center", x=0.5,title=None
        ),
        barmode='stack',
        template = "simple_white",
        xaxis=dict(tickangle=45, title=None),
        height=600,
        legend_title_text=rank.capitalize()
    )

    fig.update_xaxes(automargin=True)
    fig.update_yaxes(automargin=True)
    fig.write_html(f"{output_name}.html")
    fig.write_image(f"{output_name}.png", scale=3)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Parse final outputs to estimate abundance at a particular phylogenetic rank and generate comparative plots"
    )
    parser.add_argument("output_files", type=str,
        help="A file containing the paths of all the final species_under_genus.tsv",
    )
    parser.add_argument("rank", type=str, choices=["phylum", "class", "order", "family", "genus"],
        help="The Phylogenetic rank to group the abundances: permitted ranks are genus, family, order, class, phylum",
    )
    parser.add_argument("abundance_by", type=str,  choices=["species", "genus"], help="Genus or species level abundance")

    args = parser.parse_args()
    paths = args.output_files
    rank = args.rank
    abundance_by = args.abundance_by

    start_time = time.time() 
    modern_dfs = {}
    ancient_dfs = {}

    try:
        import kaleido
    except ImportError:
        print("You need to install 'kaleido' to export images: pip install -U kaleido")
        exit(1)

    if not os.path.exists(paths):
        print(f"Error! Please give valid file containing the output paths")
        exit(1)
    else:
        with open(paths, "r") as file:
            for path in file:
                if not os.path.exists(path.strip()):
                    print(f"WARNING: {path} doesn't exist - omitted from the analysis")
                else:
                    df = pd.read_csv(path.strip(), sep="\t")
                    df[f"{rank}"] = df["scientific_name"].apply(
                        lambda x: rank_name1(ast.literal_eval(x)[0], rank) if x else "unknown"
                        )
                    sample_name = os.path.basename(os.path.dirname(path.strip()))
                    modern_df = df[(df["score"] > 0) & (df["score"] < 4)]
                    modern_dfs[sample_name] = grouping_by_abundance(modern_df, rank, abundance_by)
                    ancient_df = df[df["score"] > 3]
                    ancient_dfs[sample_name] = grouping_by_abundance(ancient_df, rank, abundance_by)

        ancient_data = merge_dfs(ancient_dfs, rank)
        ancient_data.to_csv("ancient_abundance_matrix.tsv", sep="\t", index=False)
        heatmap_plotly(ancient_data, "ancient_taxa_abundance_heatmap", rank, abundance_by)
        barplot_plotly(ancient_data, "ancient_taxa_abundance_barplot", rank, abundance_by)

        modern_data = merge_dfs(modern_dfs, rank)
        modern_data.to_csv("modern_abundance_matrix.tsv", sep="\t", index=False)

        

    end_time = time.time()
    time_taken = end_time - start_time

    print("\n=== Abundance matrices (modern & ancient) are exported === \n")
