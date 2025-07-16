import pandas as pd
from scipy.stats import ttest_ind
import seaborn as sns
import matplotlib.pyplot as plt
import statannot


# Input data from the Pan trans paper = W10_Orthogroups.full.tsv

# Cadenza RNA seq data from Pan Tran. = 


def filter_orthogroups(file_path):
    """
    Reads a TSV file with orthogroup gene lists and filters rows where
    columns 'CS11A', 'CS11B', and 'CS11D' each contain exactly one gene (no commas).

    Parameters:
    -----------
    file_path : str
        Path to the W10_Orthogroups.full.tsv file.

    Returns:
    --------
    pandas.DataFrame
        Filtered DataFrame with columns: Orthogroup, CS11A, CS11B, CS11D
    """

    # Read the TSV file
    df = pd.read_csv(file_path, sep='\t')

    # Only keep relevant columns
    cols_to_keep = ['Orthogroup', 'CS11A', 'CS11B', 'CS11D']
    df_subset = df[cols_to_keep].copy()

    # Function to check if a cell has exactly one gene (i.e., no commas)
    def is_single_gene(cell):
        return isinstance(cell, str) and ',' not in cell and cell.strip() != ''

    # Filter rows where all three gene columns have exactly one gene
    filtered_df = df_subset[
        df_subset['CS11A'].apply(is_single_gene) &
        df_subset['CS11B'].apply(is_single_gene) &
        df_subset['CS11D'].apply(is_single_gene)
    ]

    print(f"Filtered {filtered_df.shape[0]} orthogroups with exactly one gene per CS11 genome copy.")
    return filtered_df


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

def expressed_balanced_and_stable_triads(global_30let_biases_df,  orthogroup_df, version_table, summary_df):
    """
    1. Filter global_30let_biases_df for Stable in class 
    2. Change global_30let_biases_df to be a list of the triad gene names
    3. Change gene names to v 2
    2. Filter global_30let_biases_df for their presence in summary_df in each of the columns
    """
    stable_30_lets = global_30let_biases_df[global_30let_biases_df['class'] == 'Stable'].copy()
    print(f'The stable ')
    # map the gene names to the corresponding gene names in the orthogroup_df group_id,
    # replace the column group_id with the gene names in the orthogroup_df
    # something like -stable_30_lets['group_id'] = stable_30_lets['group_id'].map(orthogroup_df.set_index('Orthogroup')['CS11A']['CS11B']['CS11D'].to_dict())
    stable_30_lets['A'] = stable_30_lets['group_id'].map(orthogroup_df.set_index('Orthogroup')['CS11A'].to_dict())
    stable_30_lets['B'] = stable_30_lets['group_id'].map(orthogroup_df.set_index('Orthogroup')['CS11B'].to_dict())
    stable_30_lets['D'] = stable_30_lets['group_id'].map(orthogroup_df.set_index('Orthogroup')['CS11D'].to_dict())
    # now convert the gene names to v2
    version_table_dict = version_table.set_index("v1.1").to_dict()["v2.1"]
    # map the gene names to the corresponding gene names in the version_table for all three columns
    stable_30_lets['A'] = stable_30_lets['A'].map(version_table_dict)
    stable_30_lets['B'] = stable_30_lets['B'].map(version_table_dict)
    stable_30_lets['D'] = stable_30_lets['D'].map(version_table_dict)
    # now filter the stable_30_lets for their presence in summary_df in each of the columns
    stable_30_lets = stable_30_lets[stable_30_lets['A'].isin(summary_df['Gene_name'])]
    stable_30_lets = stable_30_lets[stable_30_lets['B'].isin(summary_df['Gene_name'])]
    stable_30_lets = stable_30_lets[stable_30_lets['D'].isin(summary_df['Gene_name'])]
    # Get the rows where a b d are all present
    stable_30_lets = stable_30_lets[stable_30_lets['A'].notna()]
    print('stable_30_lets')
    print(stable_30_lets.head())
    pass


def read_in_expression_table_flag_leaf(expression_table_file_path):
    """
    Reads in a gene expression table from a TSV file and reformats it. 
    Parameters:
    -----------
    expression_table_file_path : str
        Path to the CS_aligned_gene_expression_table.tsv file.
    Returns:
    --------
    pandas.DataFrame
        DataFrame containing the gene expression data. e.g. - for flag leaf only 
                         gene  mean_expression_value   
                    0  TraesCS3A02G200800       1813.715101                
                    1  TraesCS5A02G358900          0.000000                
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
    expression_table = expression_table[expression_table['tissue'] == 'flag_leaf']
    # Generate mean_expression_value from the three replicates
    expression_table['mean_expression_value'] = expression_table.groupby('gene')['expression_value'].transform('mean')
    expression_table.drop(columns=['expression_value', 'tissue', 'replicate'], inplace=True)
    return expression_table

def cadenza_traids(global_30let_biases_df, orthogroup_df, version_table, expression_table, version_two_required=False):
    orthogroup_ids = global_30let_biases_df['orthogroup'].unique()
    # filter orthogroup_df for the orthogroup_idsin the orthogroup column 
    cadenza_traids = orthogroup_df[orthogroup_df['Orthogroup'].isin(orthogroup_ids)].copy()
    # Filter for  CS11A', 'CS11B', 'CS11D'
    cadenza_traids = cadenza_traids[['Orthogroup', 'CS11A', 'CS11B', 'CS11D']].copy()
    # Strip the gene name .1, .2, .3, etc. 
    if version_two_required:
        for subgenome in ['CS11A', 'CS11B', 'CS11D']:
            cadenza_traids["Gene_name"] = cadenza_traids[subgenome].str.split('.').str[0]
            cadenza_traids = convert_ref_seq_versions(cadenza_traids, version_table)
            # rename column Gene_name to the subgenome name
            # drop the subgenome column
            cadenza_traids.drop(columns=[subgenome], inplace=True)
            cadenza_traids.rename(columns={'Gene_name': subgenome}, inplace=True)

    return cadenza_traids

def assign_balance_from_expression(triads, expression_table):
    """
    A traid is unbalanced if any of the three contribute less than 20% of the total expression.
    For each row of the table, get the expression levels for the expression table 

    """
    # Strip gene names (remove .1, .2, .3 etc.)
    for col in ['CS11A', 'CS11B', 'CS11D']:
        triads[col + '_stripped'] = triads[col].str.split('.').str[0]
    gene_to_expr = dict(zip(expression_table['gene'], expression_table['mean_expression_value']))
    for col in ['CS11A', 'CS11B', 'CS11D']:
        triads[col + '_expression'] = triads[col + '_stripped'].map(gene_to_expr)

    # Calculate the total expression for each triad
    triads['total_expression'] = triads[['CS11A_expression', 'CS11B_expression', 'CS11D_expression']].sum(axis=1)
    # Calculate the percentage contribution of each gene
    triads['CS11A_percentage'] = (triads['CS11A_expression'] / triads['total_expression']) * 100
    triads['CS11B_percentage'] = (triads['CS11B_expression'] / triads['total_expression']) * 100
    triads['CS11D_percentage'] = (triads['CS11D_expression'] / triads['total_expression']) * 100

    def assign_balance(row):
        if row['CS11A_percentage'] < 20 or row['CS11B_percentage'] < 20 or row['CS11D_percentage'] < 20:
            return 'Unbalanced'
        else:
            return 'Balanced'

    triads['balance'] = triads.apply(assign_balance, axis=1)
    # Print summary statistics on number of balanced and unbalanced triads
    print("Summary of triads:")
    print(triads['balance'].value_counts())
    return triads

def assign_balance_from_exprssion_multiple_categories(triads, expression_table):
    """
    A dominant 
    A suppressed 
    B dominant
    B suppressed
    D dominant
    D suppressed
    Balanced
    """
    # Strip gene names (remove .1, .2, .3 etc.)
    for col in ['CS11A', 'CS11B', 'CS11D']:
        triads[col + '_stripped'] = triads[col].str.split('.').str[0]
    gene_to_expr = dict(zip(expression_table['gene'], expression_table['mean_expression_value']))
    for col in ['CS11A', 'CS11B', 'CS11D']:
        triads[col + '_expression'] = triads[col + '_stripped'].map(gene_to_expr)

    # Calculate the total expression for each triad
    triads['total_expression'] = triads[['CS11A_expression', 'CS11B_expression', 'CS11D_expression']].sum(axis=1)
    # Calculate the percentage contribution of each gene
    triads['CS11A_percentage'] = (triads['CS11A_expression'] / triads['total_expression']) * 100
    triads['CS11B_percentage'] = (triads['CS11B_expression'] / triads['total_expression']) * 100
    triads['CS11D_percentage'] = (triads['CS11D_expression'] / triads['total_expression']) * 100
    # if total expression is 0, skip the row
    triads = triads[triads['total_expression'] > 0]

    def assign_balance_categories(row):
        if row['CS11A_percentage'] > 60:
            return 'A dominant'
        elif row['CS11B_percentage'] > 60:
            return 'B dominant'
        elif row['CS11D_percentage'] > 60:
            return 'D dominant'
        elif row['CS11A_percentage'] > 20 and row['CS11B_percentage'] > 20 and row['CS11D_percentage'] > 20:
            return 'Balanced'
        elif row['CS11A_percentage'] > 20 and row['CS11B_percentage'] > 20:
            return 'D suppressed'
        elif row['CS11A_percentage'] > 20 and row['CS11D_percentage'] > 20:
            return 'B suppressed'
        elif row['CS11B_percentage'] > 20 and row['CS11D_percentage'] > 20:
            return 'A suppressed'
        else:
            print("Error: Triad does not fit any category")
            print(row)
            exit()
    triads['Balance_categorical'] = triads.apply(assign_balance_categories, axis=1)
    # Print summary statistics on number of balanced and unbalanced triads
    print("Summary of triads:")
    print(triads['Balance_categorical'].value_counts())
    return triads


    

def number_tcs_between_balance_unbalance(triad_table, summary_df, version_table, version_two_required=False):
    """
    For each triad, get the number of TCS between the three genes in the triad.
    """
    # Get the gene names for each triad
    triad_table['CS11A'] = triad_table['CS11A'].str.split('.').str[0]
    triad_table['CS11B'] = triad_table['CS11B'].str.split('.').str[0]
    triad_table['CS11D'] = triad_table['CS11D'].str.split('.').str[0]
    if version_two_required:
        for subgenome in ['CS11A', 'CS11B', 'CS11D']:
            triad_table["Gene_name"] = triad_table[subgenome].str.split('.').str[0]
            triad_table = convert_ref_seq_versions(triad_table, version_table)
            # rename column Gene_name to the subgenome name
            # drop the subgenome column
            triad_table.drop(columns=[subgenome], inplace=True)
            triad_table.rename(columns={'Gene_name': subgenome}, inplace=True)

    for subgenome in ['CS11A', 'CS11B', 'CS11D']:
        triad_table[subgenome + 'TC'] = triad_table[subgenome].map(summary_df.set_index('Gene_name')['Number_TC'])
    
    # Remove rows where number TC is NaN
    triad_table.dropna(subset=['CS11ATC', 'CS11BTC', 'CS11DTC'], inplace=True)

    # Total TCs
    triad_table['Total_TCs'] = triad_table[['CS11ATC', 'CS11BTC', 'CS11DTC']].sum(axis=1)
    print("\nTotal TCs by balance state:")
    print(triad_table.groupby('balance')['Total_TCs'].agg(['mean', 'std']))
    balanced = triad_table[triad_table['balance'] == 'Balanced']['Total_TCs']
    unbalanced = triad_table[triad_table['balance'] == 'Unbalanced']['Total_TCs']
    stat, p = ttest_ind(balanced, unbalanced)
    print(f"T-test p-value (Total TCs): {p}")

    # TC variation
    triad_table['TC_variation'] = (
        abs(triad_table['CS11ATC'] - triad_table['CS11BTC']) +
        abs(triad_table['CS11ATC'] - triad_table['CS11DTC']) +
        abs(triad_table['CS11BTC'] - triad_table['CS11DTC'])
    )
    print("\nTC variation by balance state:")
    print(triad_table.groupby('balance')['TC_variation'].agg(['mean', 'std']))
    balanced = triad_table[triad_table['balance'] == 'Balanced']['TC_variation']
    unbalanced = triad_table[triad_table['balance'] == 'Unbalanced']['TC_variation']
    stat, p = ttest_ind(balanced, unbalanced)
    print(f"T-test p-value (TC variation): {p}")


    def plot_distribution(data, y, y_label, title, filename):
        plt.figure(figsize=(5.7, 5.7), dpi=900)
        ax = sns.barplot(
            data=data, y=y, x='balance', hue='balance',
            errorbar='ci', capsize=0.1, palette=['#5688F7', '#E65A5C'], 
            dodge = False
        )
        # Annoate statistical significance
        statannot.add_stat_annotation( ax,
            data=data, x='balance', y=y,
            box_pairs=[("Balanced", "Unbalanced")],
            test='t-test_ind', text_format='star', loc='outside'
        )
        plt.legend().remove()
        plt.title(title, fontsize=11)
        plt.xlabel(y.replace('_', ' '), fontsize=11)
        plt.ylabel(y_label, fontsize=11)
        sns.despine()
        plt.savefig(filename, dpi=900)

    plot_distribution(
        triad_table, y='Total_TCs', y_label="Number of TSS per triad",
        title='',
        filename='total_tcs_distribution.png')
    
    plot_distribution(
        triad_table, y='TC_variation', y_label ='Variation in TSS per triad (absolute TSS count difference)',
        title='',
        filename='tc_variation_distribution.png'
    )




if __name__ == "__main__":
    orthogroup_table = "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/W10_Orthogroups.full.tsv"
    orthogroup_df = pd.read_csv(orthogroup_table, sep='\t')
    print(f'orthogroup_df = {orthogroup_df.head()}')
    corresponding_v1_v2 = "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/iwgsc_refseq_all_correspondances.csv"
    corresponding_df = pd.read_csv(corresponding_v1_v2, sep=' ')
    print(f'orthogroup_df = {corresponding_df.head()}')
    global_30let_biases = "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/global_30let_biases_table.csv"    
    global_30let_biases_df = pd.read_csv(global_30let_biases)
    print(f'global_30let_biases_df = {global_30let_biases_df.head()}')
    ## Corresponding - iwgsc_refseq_all_correspondances.csv
    ## /Volumes/Projects-Scratch/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/repos/CAGE/intermediate_data/snakemake_intermediate_data/global_30let_biases_table.csv
    summary_df = pd.read_csv("../../../intermediate_data/snakemake_intermediate_data/summary_df.csv", sep=",") # This is generated from script number_tss_to_expression_level.py
    print(f'The CAGE data is in summary df = {summary_df.head()}')
    # remove rows from summary_df where Gene_name is NaN
    summary_df.dropna(subset=['Gene_name'], inplace=True)
    #expressed_balanced_and_stable_triads(global_30let_biases_df,  orthogroup_df, corresponding_df, summary_df)

    # Generate balanced/unbalanced from cadenza RNA expression
    expression_table = read_in_expression_table_flag_leaf("/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Wheat_expression_prediction/input_data/CS_aligned_gene_expression_table.tsv")
    print(f'The expression table is = {expression_table.head()}')
    cadenza_traids = cadenza_traids(global_30let_biases_df, orthogroup_df, corresponding_df, expression_table)
    print(f'The cadenza triads are = {cadenza_traids.head()}')
    print(cadenza_traids.head())
    #cadenza_traids = assign_balance_from_expression(cadenza_traids, expression_table)
    #cadenza_traids.to_csv("/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/cadenza_triads_from_expression_rnaseq.csv", index=False)

    cadenza_triads_categorical = assign_balance_from_exprssion_multiple_categories(cadenza_traids, expression_table)
    print(f'The cadenza triads are = {cadenza_triads_categorical.head()}')
    cadenza_triads_categorical.to_csv("/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/cadenza_triads_from_expression_rnaseq_categorical.csv", index=False)

    #cadenza_triads = pd.read_csv("/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/cadenza_triads_from_expression_rnaseq.csv")
    #print(f'The cadenza triads are = {cadenza_triads.head()}')
    #number_tc = number_tcs_between_balance_unbalance(cadenza_triads, summary_df, corresponding_df, version_two_required=True)
