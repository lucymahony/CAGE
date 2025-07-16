# Script used to answer the following questions:

# Number of TC's per gene
# The relationship between the number of TSSs and the expression level of a gene

# Note: This script relies on the input CSV file being the final TC csv file - Need to think about exactly how to generate/filter this input file.

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from io import StringIO
import numpy as np
from matplotlib.ticker import ScalarFormatter
import sys





def process_cage_csv(file_path):
    """
    Processes a CAGE CSV file to summarize expression at promoter regions by gene.
    
    Parameters:
    -----------
    file_path : str
        Path to the CAGE annotation CSV file (e.g., "no_merge_distclu_md20_single_5_TC_GR_CL1_my_annot.csv").

    Returns:
    --------
    pandas.DataFrame
        A DataFrame with the following columns:
        - Gene_name: unique gene identifiers (without version suffixes)
        - Number_TC: number of tag clusters (TCs) annotated as 'promoter' for each gene
        - Expression_score: total expression score summed across those promoter TCs
    """
    df = pd.read_csv(file_path, index_col=False)
    print(f"Loaded CSV with shape: {df.shape}")
    print(f"Filtered promoter DataFrame shape: {df.shape}")
    print(df.head(), "\n")
    print(df.columns)
    df_promoter = df[df['annotation'] == 'promoter'].copy()
    # Filter to remove TC under 0.5 
    df_promoter = df_promoter[df_promoter['score'] >= 0.5]



    assert df_promoter['annotation'].nunique() == 1
    df_promoter['gene_clean'] = df_promoter['nearest_gene'].astype(str).str.split('.').str[0]

    # Group and aggregate
    summary_df = (
        df_promoter.groupby('gene_clean')
        .agg(
            Number_TC=('gene_clean', 'count'),
            Expression_score=('score', 'sum')
        )
        .reset_index()
        .rename(columns={'gene_clean': 'Gene_name'})
    )

    print("Summary DataFrame:")
    print(summary_df.head(), "\n")
    print(f"Summary DataFrame shape: {summary_df.shape}")
    assert summary_df['Gene_name'].nunique() == summary_df.shape[0], " Duplicate gene names found in summary DataFrame."
    return summary_df


def plot_number_tcs(summary_df, output_file_path):
    """
    ========== Histogram of Number_TC per Gene ==========
    """
    fig, ax = plt.subplots(figsize=(5.7, 5.7), dpi=900)
    bins = np.arange(0.5, 10.6, 1)
    sns.histplot(
        summary_df['Number_TC'],
        bins=bins,
        shrink=0.8,
        color='#0f62f3',
        edgecolor='none')
    ax.set_xlim(0.5, 10.5)
    ax.set_xticks(range(1, 11))
    # make the y-axis log
    ax.set_yscale('log')
    ax.yaxis.set_major_formatter(ScalarFormatter()) # remove scientific notation
    ax.set_xlabel("Number of Tag Clusters per Gene", fontsize=11)
    ax.set_ylabel("Number\nof\nGenes", fontsize=11, rotation=0)
    # Rotate the y axis labels
    sns.despine()
    plt.tight_layout()
    plt.savefig(output_file_path, dpi=900)

def plot_number_tss_expression_level(summary_df, output_file_path):
    """
    Boxplot of Expression Score by Number of TCs with optional trend overlay
    """
    fig, ax = plt.subplots(figsize=(5.7, 5.7), dpi=900)
    sns.set(style="whitegrid")
    palette_60 = ['#da1e28', '#d02670', '#8a3ffc', '#0f62fe', '#0072c3', '#007d79', '#198038', '#4d5358', "#726e6e"]  # Dark
    sns.boxplot(
        x='Number_TC',
        y='Expression_score',
        data=summary_df,
        ax=ax,
        showfliers=False,  # optionally hide outliers
        palette=palette_60
    )
    ax.set_title("", fontsize=11)
    ax.set_xlabel("Number of TSS", fontsize=11)
    ax.set_ylabel("Expression Score (TPM)")
    ax.set_yscale('log')
    plt.tight_layout()
    sns.despine()
    plt.savefig(output_file_path, dpi=900)



if __name__ == "__main__":
    my_annotation_file = str(sys.argv[1]) 
    output_csv_file = str(sys.argv[2])
    output_directory = str(sys.argv[3])

    summary_df = process_cage_csv(my_annotation_file)
    summary_df.to_csv(output_csv_file, index=False)
    df = pd.read_csv(output_csv_file)
    df = df.dropna() 
    top_genes = df.sort_values(by='Number_TC', ascending=False)
    top_gene = top_genes.iloc[0]
    print(f"Top gene: {top_gene['Gene_name']} ({top_gene['Number_TC']} TCs)")
    second_gene = top_genes.iloc[1]
    print(f"Second gene: {second_gene['Gene_name']} ({second_gene['Number_TC']} TCs)")
    third_gene = top_genes.iloc[2]
    print(f"Third gene: {third_gene['Gene_name']} ({third_gene['Number_TC']} TCs)")
    fourth_gene = top_genes.iloc[3]
    print(f"Fourth gene: {fourth_gene['Gene_name']} ({fourth_gene['Number_TC']} TCs)")
    fifth_gene = top_genes.iloc[4]
    print(f"Fifth gene: {fifth_gene['Gene_name']} ({fifth_gene['Number_TC']} TCs)")
    sixth_gene = top_genes.iloc[5]
    print(f"Sixth gene: {sixth_gene['Gene_name']} ({sixth_gene['Number_TC']} TCs)")
    print(df['Number_TC'].describe())
    plot_number_tcs(df, f'{output_directory}/number_tcs_per_gene.png')
    plot_number_tss_expression_level(df, f'{output_directory}/number_tss_vs_expression_level.png')