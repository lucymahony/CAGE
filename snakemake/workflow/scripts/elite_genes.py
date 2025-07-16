# Investigates if any of the genes from Hammond Kosack are in my CAGE dataset. 
import pandas as pd
import numpy as np

def read_in_version_table(version_table_file_path):
    version_table = pd.read_csv(version_table_file_path, sep=" ")
    cols_to_keep = ["v1.1", "v2.1"] 
    version_table = version_table[cols_to_keep].copy()
    return version_table

def elite_genes_table_v21(elite_genes, version_table):
    """
    Returns a df with the version 2.1 gene names 
    """
    version_table_dict = version_table.set_index("v1.1").to_dict()["v2.1"]
    elite_genes["gene_v21"] = elite_genes["gene_ID"].map(version_table_dict)
    return elite_genes

def find_elite_genes_in_summary_df(elite_genes_v21, summary_df):
    """
    Returns a df with the elite genes that are in the summary df with the summary df info 
    """
    elite_genes_v21 = elite_genes_v21.rename(columns={"gene_v21": "Gene_name"})
    elite_genes_in_summary_df = elite_genes_v21.merge(summary_df, left_on="Gene_name", right_on="Gene_name", how="left")
    print("Elite genes in summary df size:", elite_genes_in_summary_df.shape)
    # Remove rows with NaNs
    elite_genes_in_summary_df = elite_genes_in_summary_df.dropna()
    print("Elite genes in summary df size after dropping NaNs:", elite_genes_in_summary_df.shape)
    print(elite_genes_in_summary_df)
    print(f'Out of the {elite_genes_v21.shape[0]} elite genes, {elite_genes_in_summary_df.shape[0]} are in the CAGE dataset') # Out of the 82 elite genes, 27 are in the CAGE dataset
    return elite_genes_in_summary_df


if __name__ == "__main__":
    elite_genes_table = pd.read_csv("../../../intermediate_data/snakemake_intermediate_data/supplemental_reformatted_Hammond_Kosack.csv")

    summary_df = pd.read_csv("../../../intermediate_data/snakemake_intermediate_data/summary_df.csv", sep=",") # This is generated from script number_tss_to_expression_level.py

    version_table_file_path = "../../../intermediate_data/snakemake_intermediate_data/iwgsc_refseq_all_correspondances.csv"
    version_table = read_in_version_table(version_table_file_path)
    elite_genes_v21 = elite_genes_table_v21(elite_genes_table, version_table)
    elite_genes_with_cage_data = find_elite_genes_in_summary_df(elite_genes_v21, summary_df)
    elite_genes_with_cage_data.to_csv("../../../intermediate_data/snakemake_intermediate_data/elite_genes_with_cage_data.csv", index=False)

    # Are they balanced or unbalanced? 
    balance = '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/cadenza_triads_from_expression_rnaseq_categorical.csv'
    #grep "TraesCS2A02G295400" ../../../intermediate_data/snakemake_intermediate_data/cadenza_triads_from_expression_rnaseq_categorical.csv
    #OG0010903,TraesCS2A02G295400.1,TraesCS2B02G311900.1,TraesCS2D02G293200.1,TraesCS2A02G295400,TraesCS2B02G311900,TraesCS2D02G293200,35.38809672581623,17.3360602111844,0.0,52.72415693700063,67.11932211282388,32.88067788717612,0.0,A dominant
    #(base) [mahony@EI-HPC interactive scripts]$ grep "TraesCS1A02G088000" ../../../intermediate_data/snakemake_intermediate_data/cadenza_triads_from_expression_rnaseq_categorical.csv
    #OG0008208,TraesCS1A02G088000.2,TraesCS1B02G107000.3,TraesCS1D02G089300.2,TraesCS1A02G088000,TraesCS1B02G107000,TraesCS1D02G089300,1699.0819085107632,2155.4373884943902,2213.8493491718095,6068.368646176963,27.998989639187048,35.51922294391761,36.481787416895344,Balanced




