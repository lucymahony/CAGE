import pandas as pd
from typing import List

def load_data(filepath: str) -> pd.DataFrame:
    df = pd.read_csv(filepath, sep=",", index_col=False, header=None, low_memory=False, skiprows=1)
    
    columns = ['id', 'seqname', 'start', 'end', 'score', 'strand', 'dom_start', 'dom_end',
               'annotation', 'nearest_gene', 'blank', 'blank2', 'CL', 'CR', 'SIS', 'SLE', 'SRO', 'SSP']
    df.columns = columns
    print(df.columns)
    print(df.head())
    df = df[df['annotation'] == 'promoter']
    # Describe CL CR SIS SLE SRO SSP
    print(df.columns)
    print('CL', df['CL'].value_counts())
    print('CR', df['CR'].value_counts())
    print('SI', df['SIS'].value_counts())
    print('SL', df['SLE'].value_counts())
    print('SR', df['SRO'].value_counts())
    print('SS', df['SSP'].value_counts())
    # Convert tissue columns to numeric (coerce errors to NaN)
    tissue_cols = ['CL', 'CR', 'SIS', 'SLE', 'SRO', 'SSP']
    df[tissue_cols] = df[tissue_cols].apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)

    return df

def get_shared(df: pd.DataFrame, tissues: List[str]) -> pd.DataFrame:
    return df[(df[tissues] == 1).all(axis=1)]

def get_unique(df: pd.DataFrame, tissue: str, others: List[str]) -> pd.DataFrame:
    return df[(df[tissue] == 1) & (df[others].sum(axis=1) == 0)]

def save_tss_table(df: pd.DataFrame, filename: str):
    df.to_csv(filename, index=False)

def save_gene_list(df: pd.DataFrame, filename: str):
    genes = df['nearest_gene'].dropna().unique()
    genes = [g for g in genes if g and g != '']
    with open(filename, 'w') as f:
        for g in sorted(set(genes)):
            f.write(g + "\n")

def analyze_tissues(output_directory, df: pd.DataFrame, tissues: List[str], prefix: str):
    print(f"\n### TSS Sharing Analysis: {', '.join(tissues)} ###")

    # Shared
    shared = get_shared(df, tissues)
    print(f"Shared in all ({len(tissues)}): {len(shared)}")
    save_tss_table(shared, f"{output_directory}/{prefix}_TSSs_shared_all.csv")
    save_gene_list(shared, f"{output_directory}/{prefix}_genes_shared_all.txt")

    # Unique to each
    for tissue in tissues:
        others = [t for t in tissues if t != tissue]
        unique = get_unique(df, tissue, others)
        print(f"Unique to {tissue}: {len(unique)}")
        save_tss_table(unique, f"{output_directory}/{prefix}_TSSs_unique_{tissue}.csv")
        save_gene_list(unique, f"{output_directory}/{prefix}_genes_unique_{tissue}.txt")

# === MAIN USAGE ===
if __name__ == "__main__":
    # Load full overlap matrix CSV
    input_directory = "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/tmp_overlap_matrix"
    output_directory = "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data"

    df = load_data(f'{input_directory}/ALL_with_overlap_matrix.csv')
    # Example 1: CL and SLE (Cadenza vs Fielder Leaf)
    analyze_tissues(output_directory, df, ['CL', 'SLE'], prefix="CL_vs_SLE")
    analyze_tissues(output_directory, df, ['CR', 'SRO'], prefix="CR_vs_SRO")

    # Example 2: Shared and unique TSSs in 4 leaf tissues
    analyze_tissues(output_directory ,df, ['SIS', 'SLE', 'SRO', 'SSP'], prefix="SIS_SLE_SRO_SSP")
