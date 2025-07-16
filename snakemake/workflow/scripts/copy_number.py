# This script answers the following questions:
# 1. Does copy number correlate with average TSS number?
# 2. Is there a pattern in number of TCs between singletons, triads, and duplicated genes?

# Input data from PanTran paper

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import statannot


def read_in_orthogroups(orthogroup_table_file_path):
    orthogroup_table = pd.read_csv(orthogroup_table_file_path, sep="\t")
    cols_to_keep = ['Orthogroup', 'CS11A', 'CS11B', 'CS11D']
    orthogroup_table = orthogroup_table[cols_to_keep].copy()
    # Replace NaNs with empty strings
    orthogroup_table.fillna('', inplace=True)
    return orthogroup_table

def read_in_version_table(version_table_file_path):
    version_table = pd.read_csv(version_table_file_path, sep=" ")
    print(version_table.head())
    cols_to_keep = ["v1.1", "v2.1"] 
    version_table = version_table[cols_to_keep].copy()
    return version_table

def convert_ref_seq_versions(gene_table, version_table):
    # Convert gene names to v1.1 and v2.1 in the gene_name, find the match in column "v1.1" and replace with "v2.1" from the version_table
    version_table_dict = version_table.set_index("v1.1").to_dict()["v2.1"]
    gene_table["Gene_name"] = gene_table["Gene_name"].map(version_table_dict)
    return gene_table


def copy_number_table(orthogroup_table, version_table):
    # Generates a table of gene_name, copy number from the orthogroup table, then updates the version of the gene names 
    output_rows = []
    for _, row in orthogroup_table.iterrows():
        # Combine all gene names into one list across CS11A, CS11B, and CS11D
        all_genes = []
        for col in row.index[1:]:  # skip Orthogroup column
            if row[col]:
                genes = [g.strip() for g in row[col].split(',')]
                all_genes.extend(genes)
        copy_number = len(all_genes) # All genes in the same orthogroup have the same copy number
        for gene in all_genes:
            output_rows.append({'Gene_name': gene, 'Copy_number': copy_number})

    gene_table = pd.DataFrame(output_rows)
    print(f'gene_table before v2 conversion = {gene_table.head()}')
    print(f'And has dimensions {gene_table.shape}')
    # Strip the gene name .1, .2, .3, etc. 
    gene_table["Gene_name"] = gene_table["Gene_name"].str.split(".").str[0]
    # Convert gene names to v1.1 and v2.1
    gene_table_v2 = convert_ref_seq_versions(gene_table, version_table)
    print(f'gene_table after v2 conversion = {gene_table_v2.head()}')
    print(f'And has dimensions {gene_table_v2.shape}')
    return gene_table_v2

def add_number_tc_to_gene_table(gene_table, summary_df):
    # Add the number of TSS to the gene table
    gene_table = gene_table.merge(summary_df, left_on="Gene_name", right_on="Gene_name", how="left")
    print("Merged table size:", gene_table.shape)
    # Remove rows with NaNs
    gene_table = gene_table.dropna()
    # print shape of gene_table
    print('After dropping NaNs, gene_table has dimensions:', gene_table.shape)
    return gene_table


def plot_copy_number_vs_tss_old(gene_table, output_file_path):
    fig, ax = plt.subplots(figsize=(5.7, 5.7), dpi=900)
    sns.boxplot(data=gene_table, x="Copy_number", y="Number_TC", ax=ax, palette="pastel")
    plt.tight_layout()
    unique_copy_numbers = sorted(gene_table["Copy_number"].unique())
    print(f'Unique copy numbers are: {unique_copy_numbers}')
    #plt.xticks(range(len(unique_copy_numbers)), unique_copy_numbers)
    plt.savefig(output_file_path, dpi=900)

def plot_copy_number_vs_tss(gene_table, output_file_path):

    bins = [0, 1, 2, 3, 5, 10, 20, 50, gene_table["Copy_number"].max() + 1]
    labels = ["1", "2", "3", "4–5", "6–10", "11–20", "21–50", "51+"]

    gene_table["Copy_number_bin"] = pd.cut(gene_table["Copy_number"], bins=bins, labels=labels, right=True)
    bin_counts = gene_table["Copy_number_bin"].value_counts().sort_index()
    fig, ax = plt.subplots(figsize=(5.7, 5.7), dpi=900)
    sns.barplot(data=gene_table,
        x="Copy_number_bin",
        y="Number_TC",
        ax=ax,
        ci="sd",
        palette="pastel")
    ax.set_xlabel("Copy number bin")
    ax.set_ylabel("Number of TCs")
    # Add "n = X" under each x-axis tick
    xtick_labels = [f"{label}\nn={bin_counts[label]}" for label in labels]
    ax.set_xticklabels(xtick_labels)
    plt.tight_layout()
    print("Genes per bin:")
    print(bin_counts)
    plt.savefig(output_file_path, dpi=900)



# Also singletons vs triads:
# Take the orthogroup_table_file_path and generate a df of single/triad genes
def singleton_vs_triad(orthogroup_table, version_table):
    output_rows = []
    # Replace NaNs with empty strings for safety
    orthogroup_table = orthogroup_table.fillna('')
    for _, row in orthogroup_table.iterrows():
        # Get gene lists from each subgenome
        a = [g.strip() for g in row["CS11A"].split(',') if g.strip()]
        b = [g.strip() for g in row["CS11B"].split(',') if g.strip()]
        d = [g.strip() for g in row["CS11D"].split(',') if g.strip()]
        total_genes = a + b + d
        num_genes = len(total_genes)
        if len(a) == 1 and len(b) == 1 and len(d) == 1:
            status = "triad"
        elif num_genes == 1:
            status = "singleton"
        else:
            status = "duplicated"
        
        for gene in total_genes:
            output_rows.append({'Gene_name': gene, 'Copy_status': status})
    gene_table = pd.DataFrame(output_rows)
    gene_table["Gene_name"] = gene_table["Gene_name"].str.split('.').str[0]
    gene_table_v2 = convert_ref_seq_versions(gene_table, version_table)
    return gene_table_v2

def plot_triad_cn(gene_table, output_file_path):
    # print value counts for gene_table
    print(gene_table["Copy_status"].value_counts())
    fig, ax = plt.subplots(figsize=(5.7, 5.7), dpi=900)
    order = ["singleton", "triad", "duplicated"]
    palette_60 = ['#da1e28', '#d02670', '#8a3ffc']
    sns.barplot(data=gene_table, x="Copy_status", y="Number_TC", ax=ax, errorbar="ci", capsize=0.1, palette=palette_60, order=order)
    statannot.add_stat_annotation( ax,
            data=gene_table, x='Copy_status', y="Number_TC",
            box_pairs=[("singleton", "triad"), ("singleton", "duplicated"),("triad", "duplicated") ],
            test='Kruskal', text_format='star', loc='outside'
        )
    ax.set_xlabel("Gene copy number status")
    ax.set_ylabel("Number of TSS per gene")
    plt.tight_layout()
    # despine the plot
    sns.despine()
    plt.savefig(output_file_path, dpi=900)

if __name__ == "__main__":
    orthogroup_table_file_path = "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/W10_Orthogroups.full.tsv"
    version_table_file_path = "../../../intermediate_data/snakemake_intermediate_data/iwgsc_refseq_all_correspondances.csv"
    summary_df = pd.read_csv("../../../intermediate_data/snakemake_intermediate_data/summary_df.csv", sep=",") # This is generated from script number_tss_to_expression_level.py
    orthogroup_table = read_in_orthogroups(orthogroup_table_file_path)
    version_table = read_in_version_table(version_table_file_path)
    gene_table_v2_cn = copy_number_table(orthogroup_table, version_table)
    gene_table_v2_cn = add_number_tc_to_gene_table(gene_table_v2_cn, summary_df)
    plot_copy_number_vs_tss(gene_table_v2_cn, "../../../intermediate_data/snakemake_intermediate_data/copy_number_vs_tss.png")
    triad_cn = singleton_vs_triad(orthogroup_table, version_table)
    triad_cn = add_number_tc_to_gene_table(triad_cn, summary_df)
    plot_triad_cn(triad_cn, "../../../intermediate_data/snakemake_intermediate_data/triad_vs_tss.png")