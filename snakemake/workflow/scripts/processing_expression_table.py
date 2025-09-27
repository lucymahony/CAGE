# Script to take CS_aligned_gene_expression_table.tsv and generate a dataset of gene_name, expression_value, tissue 

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import norm
from sklearn.cluster import KMeans

def read_in_expression_table(expression_table_file_path):
    """
    Reads in a gene expression table from a TSV file and reformats it. 
    Parameters:
    -----------
    expression_table_file_path : str
        Path to the CS_aligned_gene_expression_table.tsv file.
    Returns:
    --------
    pandas.DataFrame
        DataFrame containing the gene expression data. e.g. 
                         gene  expression_value tissue  replicate
                    0  TraesCS3A02G200800       1813.715101      D          1
                    1  TraesCS5A02G358900          0.000000      D          1
    """
    expression_table = pd.read_csv(expression_table_file_path, sep="\t")
    #D: dawn, F: flag leaf, G: grain, R: root, S: spike, V: dusk.
    tissue_column_names = [f"CAD_{tissue}{replicate}" for tissue in ['D', 'F', 'G', 'R', 'S', 'V'] for replicate in range(1, 4)]
    cols_to_keep = ['gene'] + tissue_column_names
    expression_table = expression_table[cols_to_keep].copy()
    expression_table = expression_table.melt(id_vars=['gene'], var_name='tissue_replicate', value_name='expression_value')
    expression_table[['tissue', 'replicate']] = expression_table['tissue_replicate'].str.extract(r'CAD_(\w)(\d)')
    expression_table.drop(columns=['tissue_replicate'], inplace=True)
    expression_table['replicate'] = expression_table['replicate'].astype(int)
    expression_table['expression_value'] = expression_table['expression_value'].astype(float)
    expression_table.dropna(subset=['expression_value'], inplace=True)
    # Replace codes with full names
    tissue_mapping = {
        'D': 'dawn',
        'F': 'flag_leaf',
        'G': 'grain',
        'R': 'root',
        'S': 'spike',
        'V': 'dusk'
    }
    expression_table['tissue'] = expression_table['tissue'].map(tissue_mapping)
    print(expression_table.head())
    # transform the expression values to log2 +1 # new column 
    expression_table['log_expression_value'] = expression_table['expression_value'].apply(lambda x: np.log2(x + 1))
    return expression_table


def plot_replicates(expression_table, tissue, rep_number_1, rep_number_2, output_file_path, log=False):
    """
    Plots a scatter plot of the expression values of two biosample numbers (2 replicates) from the same study / the same expression matrix. 
    """
    fig, ax = plt.subplots(figsize=(5.7, 5.7))
    # expression_table 
    ax.scatter(expression_table[expression_table['tissue'] == tissue]['expression_value'][expression_table['replicate'] == rep_number_1],
                expression_table[expression_table['tissue'] == tissue]['expression_value'][expression_table['replicate'] == rep_number_2],
                color='blue')
    if log:
        ax.set_xscale("log")
        ax.set_yscale("log")
    # X label is tissue plus replicate number
    plt.xlabel("Replicate 2 Expression Level TPM")
    plt.ylabel("Replicate 1 Expression Level TPM")
    # Display the plot
    plt.savefig(output_file_path, dpi=900)



def ma_plot_boxplot(expression_table, out_file_path, title='MA Plot'):
    """
    Given an experiment dictionary it plots the 2 biological replicates in the style of an MA plot.
    x-axis = log average = A = ( log(a) + log(b) )/ 2
    y-axis = log ratio = M = (log(a/b)) = log(a) - log(b)
    where:
    a = rep 1
    b = rep 2
    """
    # Make a table of gene_name, log_average, log_ratio
    ma_table = expression_table[['gene', 'tissue', 'replicate', 'log_expression_value']].copy()
    ma_table = ma_table[ma_table['tissue'] == 'flag_leaf']
    ma_table = ma_table.pivot(index='gene', columns='replicate', values='log_expression_value')
    ma_table = ma_table.rename(columns={1: 'log_expression_value_rep1', 2: 'log_expression_value_rep2'})
    ma_table = ma_table.reset_index()
    print('MA table:')
    print(ma_table.head())

    ma_table['log_average'] = (ma_table['log_expression_value_rep1'] + ma_table['log_expression_value_rep2']) / 2
    ma_table['log_ratio'] = ma_table['log_expression_value_rep1'] - ma_table['log_expression_value_rep2']
    # Create bins of size 1 for A
    ma_table['A_bin'] = np.floor(ma_table['log_average']).astype(int)

    # Create the plot with box plots for each bin of A
    fig, ax = plt.subplots(figsize=(5.7, 5.7), dpi=900)
    sns.boxplot(x='A_bin', y='log_ratio', data=ma_table, ax=ax, color='#0f62fe')
    ax.axhline(y=1, color='grey', linestyle='--', zorder=10)
    ax.axhline(y=-1, color='grey', linestyle='--', zorder=20)
    ax.set_xlabel('A (Log2 Average)', fontsize=11)
    ax.set_ylabel('M (Log2 Ratio)', fontsize=11)
    ax.set_yticks(np.linspace(-6, 6, num=7))   
    ax.set_title(title, fontsize=11)
    fig.savefig(out_file_path, dpi=900)

def generate_average_log_expression_table(expression_table_file_path, output_file_path):
    expression_table = pd.read_csv(expression_table_file_path, sep="\t")
    #D: dawn, F: flag leaf, G: grain, R: root, S: spike, V: dusk.
    tissue_column_names = [f"CAD_{tissue}{replicate}" for tissue in ['D', 'F', 'G', 'R', 'S', 'V'] for replicate in range(1, 4)]
    cols_to_keep = ['gene'] + tissue_column_names
    expression_table = expression_table[cols_to_keep].copy()
    expression_table = expression_table.melt(id_vars=['gene'], var_name='tissue_replicate', value_name='expression_value')
    expression_table[['tissue', 'replicate']] = expression_table['tissue_replicate'].str.extract(r'CAD_(\w)(\d)')
    expression_table.drop(columns=['tissue_replicate'], inplace=True)
    expression_table['replicate'] = expression_table['replicate'].astype(int)
    expression_table['expression_value'] = expression_table['expression_value'].astype(float)
    expression_table.dropna(subset=['expression_value'], inplace=True)
    # Replace codes with full names
    tissue_mapping = {
        'D': 'dawn',
        'F': 'flag_leaf',
        'G': 'grain',
        'R': 'root',
        'S': 'spike',
        'V': 'dusk'
    }
    expression_table['tissue'] = expression_table['tissue'].map(tissue_mapping)
    expression_table = expression_table[expression_table['tissue'] == 'flag_leaf']
    # df of gene_name, average_log_expression_value
    expression_table = expression_table.groupby('gene')['expression_value'].mean().reset_index()
    expression_table = expression_table[expression_table['expression_value'] >= 0.5]
    expression_table['log_expression_value'] = expression_table['expression_value'].apply(lambda x: np.log2(x + 1))
    expression_table.to_csv(output_file_path, sep=",", index=False)
    print(f'Average log2 x+1 expression table written to: {output_file_path}')



def compute_tau(expression_table, tissues, tau_threshold=0.85, z=1.0) -> pd.DataFrame:
    """
    Computes tissue specificity (tau) per gene across selected tissues.

    Genes with expression < 1 TPM in all tissues get tissue = 'null' and tau = NaN.
    Genes with tau < threshold are labeled as 'widespread'.
    Tau calculated by     
        1) filtering out genes with mean expression < 1 TPM in any tissue
        2) log2-transforming remaining values. log2(x+1) the x +1 required so that no negative values. 
        3) Compute tissue-specificity (tau)     

    Parameters
    ----------
    expression_table : pd.DataFrame
        Columns: 'gene', 'expression_value', 'tissue', 'replicate'
    tissues : list
        List of tissue names to include.
    tau_threshold : float
        Tau threshold for tissue-specificity classification.
    z : float
        Z-score cutoff for calling tissue specificity.

    Returns
    -------
    pd.DataFrame with columns: 'gene', 'tissue', 'tau'
    """
    expression_table = expression_table[expression_table['tissue'].isin(tissues)]
    mean_expr = expression_table.groupby(['gene', 'tissue'])['expression_value'].mean().reset_index()

    expr_wide = mean_expr.pivot(index='gene', columns='tissue', values='expression_value')
    expr_wide = expr_wide.dropna()

    expr_wide_filtered = expr_wide[(expr_wide >= 1).any(axis=1)]
    filtered_out_genes = expr_wide[~(expr_wide >= 1).any(axis=1)]


    def calculate_tau(x):
        xmax = x.max()
        if xmax == 0 or np.isnan(xmax):
            return np.nan
        x_hat = x / xmax
        tau = ((1 - x_hat).sum()) / (len(x_hat) - 1)
        return tau

    expr_log2 = np.log2(expr_wide_filtered + 1)
    tau_values = expr_log2.apply(calculate_tau, axis=1)
    expr_log2['tau'] = tau_values

    # 1. Tissue-specific genes (tau >= threshold)
    specific_candidates = expr_log2[expr_log2['tau'] >= tau_threshold]
    print(f"{specific_candidates.shape[0]} genes pass τ ≥ {tau_threshold} filter")

    tissue_specific_rows = []
    for gene, row in specific_candidates.drop(columns='tau').iterrows():
        xmax = row.max()
        sigma = row.std()
        if pd.isna(sigma) or sigma == 0:
            continue

        specific_tissues = [tissue for tissue, value in row.items() if value >= xmax - (z * sigma)]


        if len(specific_tissues) == 0 or len(specific_tissues) == len(row):
            continue

        for tissue in specific_tissues:
            tissue_specific_rows.append({
                'gene': gene,
                'tissue': tissue,
                'tau': specific_candidates.loc[gene, 'tau']
            })

    result_df = pd.DataFrame(tissue_specific_rows)

    # 2. Widespread genes (tau < threshold)
    widespread_genes = expr_log2[expr_log2['tau'] < tau_threshold].copy()
    widespread_df = widespread_genes[['tau']].reset_index()
    widespread_df['tissue'] = 'widespread'

    # 3. Null genes (filtered out for low expression)
    null_df = pd.DataFrame({
        'gene': filtered_out_genes.index,
        'tissue': 'null',
        'tau': np.nan
    })

    # Combine all
    full_result = pd.concat([result_df, widespread_df, null_df], ignore_index=True)
    n_null = (full_result['tissue'] == 'null').sum()
    n_widespread = (full_result['tissue'] == 'widespread').sum()
    n_specific = (full_result['tissue'] != 'null') & (full_result['tissue'] != 'widespread')
    specific_counts = full_result[n_specific].groupby('gene').size()
    n_specific = specific_counts.count()
    n_specific_single = (specific_counts == 1).sum()

    print(f'There are {n_null} genes assigned as null.')
    print(f'There are {n_widespread} genes assigned as widespread.')
    print(f'There are {n_specific} genes assigned as tissue-specific.')
    print(f'There are {n_specific_single} genes assigned as tissue-specific in only one tissue.')

    return mean_expr, full_result



def compute_cluster_centers(expression_matrix: pd.DataFrame) -> pd.DataFrame:
    centers = []
    for gene, row in expression_matrix.iterrows():
        values = row.values.reshape(-1, 1)
        if np.count_nonzero(values) < 2:
            centers.append([np.nan, np.nan])
            continue
        kmeans = KMeans(n_clusters=2, random_state=0, n_init='auto').fit(values)
        centers.append(sorted(kmeans.cluster_centers_.flatten()))
    return pd.DataFrame(centers, columns=['lower_cluster', 'upper_cluster'], index=expression_matrix.index)

def compute_threshold_ratios(expression_matrix: pd.DataFrame, cluster_centers: pd.DataFrame) -> pd.Series:
    ratios = []
    for gene in expression_matrix.index:
        expr_values = expression_matrix.loc[gene].values
        nonzero = expr_values[expr_values > 0]
        if len(nonzero) == 0 or np.isnan(cluster_centers.loc[gene, 'lower_cluster']):
            ratios.append(np.nan)
            continue
        cluster_mean = cluster_centers.loc[gene].mean()
        above_mean = np.sum(expr_values > cluster_mean)
        ratios.append(above_mean / len(nonzero))
    return pd.Series(ratios, index=expression_matrix.index, name='threshold_ratio')

def compute_z_threshold_from_ratios(ratios: pd.Series) -> float:
    ratios_clean = ratios.dropna().reset_index(drop=True)
    x = np.arange(1, len(ratios_clean) + 1)
    coeffs = np.polyfit(x, ratios_clean, deg=1)
    optimized_ratio = np.mean(coeffs[0] * x + coeffs[1])
    z_threshold = norm.ppf(0.95) - abs(norm.ppf(optimized_ratio))
    return z_threshold

if __name__ == "__main__":
    expression_table_file_path = "../../../input_data/CS_aligned_gene_expression_table.tsv"
    expression_table = read_in_expression_table(expression_table_file_path)
    print(expression_table.head())
    print(f"The tissues are {expression_table['tissue'].unique()}")
    print(f"The tissues are {expression_table['replicate'].unique()}")
    pivoted_expression = (expression_table.groupby(['gene', 'tissue'])['log_expression_value'].mean().unstack(fill_value=0))

    #centers = compute_cluster_centers(pivoted_expression)
    #ratios = compute_threshold_ratios(pivoted_expression, centers)
    #z_val = compute_z_threshold_from_ratios(ratios) 
    #print(f' the z val is {z_val}')
    # z val is 1.468448433935507
    z_val=1.468448433935507

    mean_expr, result_df = compute_tau(expression_table, ['flag_leaf', 'grain', 'root', 'spike'], tau_threshold=0.85, z=z_val)

    # Save the average expression per gene
    # '../../../tmp/snakemake_intermediate_data'
    mean_expr.to_csv(f"../../../tmp/snakemake_intermediate_data/cadenza_mean_expression_per_tissue.csv", index=False)
     
    # Save the tissue specific genes
    result_df.to_csv(f"../../../tmp/snakemake_intermediate_data/cadenza_tissue_specific_list_z{z_val}.csv", index=False)
