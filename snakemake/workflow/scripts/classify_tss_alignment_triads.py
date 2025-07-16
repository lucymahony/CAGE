# Script to classify the alignment of TSS between genes within a triad. E.g. are they in the same place relative to the START codon 

# A) Takes at categorical approach defining triads as:
# 1) Same - all 3 in the same place
# 2) Shiffted  - two in the same place, one shifted
# 3) Different - all 3 in different places

# B) A continuous approach where for each pairwise comparison, the number of same / total pairs so 1 = perfect match and 0 = no match

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
from scipy.stats import chi2_contingency


def build_dict_start_codons(genome_annotation_file_path):
    """
    Build a dictionary of start codons from a GFF/GTF file for each gene.
    Dictionary keys are gene IDs and values are (start_codon, strand).
    Only considers first CDS feature seen per gene.
    """
    start_codon_dict = {}

    with open(genome_annotation_file_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue  # skip headers and empty lines

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue  # malformed line, skip

            seqname, source, feature_type, start, end, score, strand, frame, attributes = fields

            if feature_type != 'CDS':
                continue  # only interested in CDS features

            # Extract gene ID from attributes
            gene_id = None
            if 'Parent=' in attributes:
                gene_id = attributes.split('Parent=')[1].split(';')[0]
            elif 'ID=' in attributes:  # fallback
                gene_id = attributes.split('ID=')[1].split(';')[0]

            if gene_id is None:
                continue  # no gene ID found, skip

            # Clean up the gene_id to match your format if needed
            gene_id = gene_id.split('.')[0]  # remove transcript/isoform part if needed

            # Only store the first CDS seen for a gene
            if gene_id not in start_codon_dict:
                start = int(start)
                end = int(end)
                if strand == '+':
                    start_codon = start
                else:
                    start_codon = end
                start_codon_dict[gene_id] = (start_codon, strand)

    return start_codon_dict

def build_dict_tss(tss_file_path):
    """
    tss_file_path looks like: 
    "seqnames","start","end","width","strand","score","dominant_ctss.seqnames","dominant_ctss.pos","dominant_ctss.strand","dominant_ctss.score","nr_ctss","annotation","genes","nearest_gene","gene"
    "Chr1A",213,214,2,"+",0.213433390906991,"Chr1A",213,"+",0.106716695453495,2,"unknown","","","1"
    "Chr1A",478,494,17,"+",1.28060034544194,"Chr1A",478,"+",0.426866781813981,5,"unknown","","","1"
    # Tss dict of gene id, tss position, strand 
    """
    # Read tss into a pandas dataframe
    tss_df = pd.read_csv(tss_file_path, sep=",", header=0)
    # filter for pomoter in "annotation" column
    tss_df = tss_df[tss_df['annotation'].str.contains("promoter")]
    print(f'There are {tss_df.shape[0]} TSS in the tss_df after filtering for promoters.')
    # Create a dictionary of TSS positions for each gene
    # For muliple tss make a list of tss positions and strand
    tss_dict = {}
    for index, row in tss_df.iterrows():
        gene_id = row['nearest_gene']
        if pd.isna(gene_id):
            continue
        gene_id = gene_id.split('.')[0]  # remove transcript/isoform part if needed
        tss_position = int(row['start'])
        strand = row['strand']
        if gene_id not in tss_dict:
            tss_dict[gene_id] = [(tss_position, strand)]
        tss_dict[gene_id].append((tss_position, strand))
    return tss_dict

def build_distance_to_start_codon_dict(start_codon_dict, tss_dict):
    """
    Build a dictionary of distances to start codons for each gene
    The dictionary keys are gene IDs and the values are lists of distances to start codons.
    """
    distance_to_start_codon_dict = {}
    for gene_id, tss_positions in tss_dict.items():
        if gene_id in start_codon_dict:
            start_codon_pos, start_codon_strand = start_codon_dict[gene_id]
            distances = []
            for tss_position, tss_strand in tss_positions:
                if start_codon_strand != tss_strand:
                    print(f"Warning: TSS and start codon strands do not match for gene {gene_id}.")
                else:
                    distance = abs(tss_position - start_codon_pos)
                distances.append(distance)
            distance_to_start_codon_dict[gene_id] = distances
    return distance_to_start_codon_dict

def generate_triad_table(orthogroup_file_path):
    # Takes orthogroup_table_file_path and generates a triad table
    df = pd.read_csv(orthogroup_file_path, sep='\t')
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
        df_subset['CS11D'].apply(is_single_gene)]
    print(f"There are {filtered_df.shape[0]} triads -e.g. orthogroups with exactly one gene per CS11 genome copy.")
    return filtered_df

def generate_triad_distance_to_start_codon_table(distance_to_start_codon_dict, triad_df, corresponding_dict):
    # Generates a df of 
    # A_gene_id, B_gene_id, D_gene_id, A_distance_to_start_codon, B_distance_to_start_codon, D_distance_to_start_codon
    # where distance to start codon is a list of ints. 
    rows = []

    for index, row in triad_df.iterrows():
        orthogroup_id = row['Orthogroup']
        
        a_gene = row['CS11A'].split('.')[0]
        b_gene = row['CS11B'].split('.')[0]
        d_gene = row['CS11D'].split('.')[0]
        #print(f' The genes are {a_gene}, {b_gene}, {d_gene}')
        # Convert 2G to 3G 
        a_corresponding = corresponding_dict.get(a_gene, None) 
        b_corresponding = corresponding_dict.get(b_gene, None)
        d_corresponding = corresponding_dict.get(d_gene, None)
        #print(f'The corresponding genes are {a_corresponding}, {b_corresponding}, {d_corresponding}')

        # Get distances. 
        a_distances = distance_to_start_codon_dict.get(a_corresponding, np.nan) # If there are no TSS for this gene it will return NaN
        b_distances = distance_to_start_codon_dict.get(b_corresponding, np.nan)
        d_distances = distance_to_start_codon_dict.get(d_corresponding, np.nan)

        #print(f'The distances are {a_distances}, {b_distances}, {d_distances}')

        # Save the row
        rows.append({
            'Orthogroup': orthogroup_id,
            'A_gene': a_gene,
            'B_gene': b_gene,
            'D_gene': d_gene,
            'A_distances': a_distances,
            'B_distances': b_distances,
            'D_distances': d_distances
        })
    print(f'There are {len(rows)} rows in the triad distance df.')
    # Remove rows with NaN distances
    triad_distance_df = pd.DataFrame(rows)
    triad_distance_df = triad_distance_df.dropna(axis=0)
    print(f'There are {triad_distance_df.shape[0]} rows in the triad distance df after dropping NaN distances.')
    return triad_distance_df

def coresponding_gene_names():
    corresponding_v1_v2 = "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/iwgsc_refseq_all_correspondances.csv"
    corresponding_df = pd.read_csv(corresponding_v1_v2, sep=' ') # v1.1, v2.1
    # Turn a df into a dict with v1.1 as the key and v2.1 as the value
    corresponding_dict = {}
    for index, row in corresponding_df.iterrows():
        v1 = row['v1.1']
        v2 = row['v2.1']
        if pd.isna(v1) or pd.isna(v2):
            continue
        v1 = str(v1).split('.')[0] 
        v2 = str(v2).split('.')[0] 
        corresponding_dict[v1] = v2
    return corresponding_dict


def classify_triad_tss_categorical(triad_distance_df, threshold):
    """
    Classify TSS alignment for each triad based on distances to start codon.
    Classifies as 'Same', 'Shifted', or 'Different'.

    For multiple TSS, if any of the TSS are the same then its counted as the same. 
    """
    def check_match(list1, list2, threshold):
        return any(abs(v1 - v2) <= threshold for v1 in list1 for v2 in list2)

    def classify(row):
        a_distances = row['A_distances']
        b_distances = row['B_distances']
        d_distances = row['D_distances']
        # Remove duplicates
        a_distances = list(set(a_distances))
        b_distances = list(set(b_distances))
        d_distances = list(set(d_distances))

        matches = sum([
            check_match(a_distances, b_distances, threshold),
            check_match(a_distances, d_distances, threshold),
            check_match(b_distances, d_distances, threshold)])
        return ("Same" if matches == 3 else "Shifted" if matches == 2 else "Different")

    triad_distance_df['TSS_alignment'] = triad_distance_df.apply(classify, axis=1)
    return triad_distance_df


def classify_triad_tss_continuous(triad_distance_df, threshold):
    def tss_alignment_score(tss_A, tss_B, tss_D, threshold=100):
        # Get all possible TSS pairs
        total_pairs = 0
        matched_pairs = 0

        if not tss_A or not tss_B or not tss_D:
            return np.nan 
        total_pairs += len(tss_A) * len(tss_B)
        matched_pairs += sum(1 for a in tss_A for b in tss_B if abs(a - b) <= threshold)
        total_pairs += len(tss_A) * len(tss_D)
        matched_pairs += sum(1 for a in tss_A for d in tss_D if abs(a - d) <= threshold)
        total_pairs += len(tss_B) * len(tss_D)
        matched_pairs += sum(1 for b in tss_B for d in tss_D if abs(b - d) <= threshold)

        if total_pairs == 0:
            return np.nan

        return matched_pairs / total_pairs
    def classify(row):
        a_distances = row['A_distances']
        b_distances = row['B_distances']
        d_distances = row['D_distances']
        # Remove duplicates
        a_distances = list(set(a_distances))
        b_distances = list(set(b_distances))
        d_distances = list(set(d_distances))
        score = tss_alignment_score(a_distances, b_distances, d_distances, threshold)
        return score
    triad_distance_df['TSS_alignment'] = triad_distance_df.apply(classify, axis=1)
    return triad_distance_df



def generate_distance_database(genome_annotation_file_path, orthogroup_table_file_path, tss_file_path):
    start_codon_dict = build_dict_start_codons(genome_annotation_file_path)
    #print(f'There are {len(start_codon_dict)} genes with start codons.')
    tss_dict = build_dict_tss(tss_file_path)
    #print(f'There are {len(tss_dict)} genes with TSS.')
    distance_to_start_codon_dict = build_distance_to_start_codon_dict(start_codon_dict, tss_dict)
    #print(f'There are {len(distance_to_start_codon_dict)} genes with TSS and start codons.')    
    #print(f'The average distance to start codon is {np.mean([np.mean(v) for v in distance_to_start_codon_dict.values()])}')
    # The average distance to start codon is 308.9856306239929
    corresponding_dict = coresponding_gene_names()
    triad_df = generate_triad_table(orthogroup_table_file_path)
    triad_distance_df = generate_triad_distance_to_start_codon_table(distance_to_start_codon_dict, triad_df, corresponding_dict)
    #print(triad_distance_df.head()) 
    #    There are 3861 rows in the triad distance df after dropping NaN distances.
    #   Orthogroup              A_gene  ...   B_distances         D_distances
    #34  OG0003029  TraesCS2A02G164000  ...    [178, 178]          [162, 162]
    #36  OG0003114  TraesCS1A02G440200  ...    [178, 178]          [208, 208]
    #69  OG0003929  TraesCS7A02G111400  ...  [4025, 4025]  [2477, 2477, 2208]
    #78  OG0004057  TraesCS6A02G113500  ...    [116, 116]          [104, 104]
    #90  OG0004296  TraesCS3A02G470700  ...    [129, 129]          [127, 127]
    return triad_distance_df


def classify_tss_alignment_categorical(triad_distance_df, threshold, balance_classification_file_path, output_df_file_path):
    """
    Classify TSS alignment for each triad based on distances to start codon.
    Classifies as 'Same', 'Shifted', or 'Different'.
    """
    triad_tss_classification = classify_triad_tss_categorical(triad_distance_df, threshold)
    print(triad_tss_classification['TSS_alignment'].value_counts())

    balance_classification_categorical = pd.read_csv(balance_classification_file_path, sep=',', header=0)

    balance_classification_categorical = balance_classification_categorical[['CS11A_stripped', 'CS11B_stripped', 'CS11D_stripped', 'Balance_categorical']]
    balance_classification_categorical = balance_classification_categorical.rename(columns={
        'CS11A_stripped': 'A_gene',
        'CS11B_stripped': 'B_gene',
        'CS11D_stripped': 'D_gene',
        'Balance_categorical': 'balance'
    })

    triad_tss_classification = triad_tss_classification.merge(
        balance_classification_categorical[['A_gene', 'balance']],
        on='A_gene',
        how='left')
    triad_tss_classification = triad_tss_classification.dropna(subset=['balance'])
    print(triad_tss_classification.head())
    print(f"There are now {triad_tss_classification.shape[0]} rows in the triad_tss_classification dataframe.")
    # Save the triad_tss_classification dataframe to a csv file
    triad_tss_classification.to_csv(output_df_file_path, index=False)


def classify_tss_alignment_continuous(triad_distance_df, threshold, balance_classification_file_path, output_df_file_path):
    """
    Classify TSS alignment for each triad based on distances to start codon.
    Classifies as a continuous score between 0 and 1.
    """
    triad_tss_classification = classify_triad_tss_continuous(triad_distance_df, threshold)
    balance_classification_categorical = pd.read_csv(balance_classification_file_path, sep=',', header=0)
    balance_classification_categorical = balance_classification_categorical[['CS11A_stripped', 'CS11B_stripped', 'CS11D_stripped', 'Balance_categorical']]
    balance_classification_categorical = balance_classification_categorical.rename(columns={
        'CS11A_stripped': 'A_gene',
        'CS11B_stripped': 'B_gene',
        'CS11D_stripped': 'D_gene',
        'Balance_categorical': 'balance'
    })
    triad_tss_classification = triad_tss_classification.merge(
        balance_classification_categorical[['A_gene', 'balance']],
        on='A_gene',
        how='left')
    triad_tss_classification = triad_tss_classification.dropna(subset=['balance'])
    print(triad_tss_classification.head())
    print(f"There are now {triad_tss_classification.shape[0]} rows in the triad_tss_classification dataframe.")
    # return summary of the data in the TSS_alignment column
    print(triad_tss_classification['TSS_alignment'].describe())
    return triad_tss_classification

def plot_continuous(triad_tss_classification, threshold):
    # Plot the TSS alignment scores
    plt.figure(figsize=(10, 6))
    sns.barplot(triad_tss_classification, x='balance', y='TSS_alignment')
    plt.title(f'TSS Alignment Scores by Balance Category (Threshold = {threshold})')
    plt.xlabel('Balance Category')
    plt.ylabel('TSS Alignment Score')
    plt.savefig(f'{threshold}_triad_tss_classification_continuous.png', dpi=900)


if __name__ == "__main__":
    # ---- Input files and parameters ---- 
    genome_annotation_file_path = sys.argv[1]  # Path to the genome annotation file e.g. iwgsc_refseqv2.1_annotation_200916_HC.gff3 
    orthogroup_table_file_path = sys.argv[2]  # Path to the triad file
    tss_file_path = sys.argv[3]  # Path to the TSS file
    threshold = int(sys.argv[4]) if len(sys.argv) > 3 else 10  # Threshold for TSS classification
    print(f'The input files are {genome_annotation_file_path}, {orthogroup_table_file_path}, {tss_file_path} and the threshold is {threshold}')

    # ---- Generate distance database  ---- 
    triad_distance_df = generate_distance_database(genome_annotation_file_path, orthogroup_table_file_path, tss_file_path)

    # ---- Classify TSS alignment Categorical ----
    balance_classification_file_path='/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/cadenza_triads_from_expression_rnaseq_categorical.csv'
    output_df_file_path = f'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/triad_tss_classification_{threshold}.csv'
    #classify_tss_alignment_categorical(triad_distance_df, threshold, balance_classification_file_path, output_df_file_path)
    #triad_tss_classification = pd.read_csv(f'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/triad_tss_classification_{threshold}.csv', sep=',', header=0)
    #pattern_table = pd.crosstab(triad_tss_classification['TSS_alignment'], triad_tss_classification['balance'])
    #print(pattern_table)
    #pattern_table.to_csv(f'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/pattern_table_{threshold}.csv', index=True)
    #pattern_table_swapped = pd.crosstab(
    #    triad_tss_classification['balance'],
    #    triad_tss_classification['TSS_alignment'])

    ## Normalize to proportions within each balance category
    #pattern_table_swapped_prop = pattern_table_swapped.div(pattern_table_swapped.sum(axis=1), axis=0)

    ## Plot
    #pattern_table_swapped_prop.plot(
    #    kind='bar', 
    #    stacked=False, 
    #    figsize=(10,6),
    #    title='Proportion of TSS Alignment Classes per Balance Category'
    #)
    #plt.ylabel('Proportion')
    #plt.xlabel('Balance Category')
    #plt.legend(title='TSS Alignment', bbox_to_anchor=(1.05, 1), loc='upper left')
    #plt.xticks(rotation=45, ha='right')
    #plt.tight_layout()
    #plt.savefig(f'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/{threshold}_triad_tss_classification.png', dpi=900)
    #chi2, p, dof, expected = chi2_contingency(pattern_table_swapped)
    #print(f"Chi-square statistic = {chi2:.2f}")
    #print(f"Degrees of freedom = {dof}")
    #print(f"P-value = {p:.4e}")

    # ---- Classify TSS alignment Continuous ----
    triad_tss_classification_continuous = classify_tss_alignment_continuous(triad_distance_df, threshold, balance_classification_file_path, output_df_file_path)
    # Optionally save the continuous classification as a csv. 
    #triad_tss_classification_continuous = pd.read_csv(f'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/triad_tss_classification_{threshold}.csv', sep=',', header=0)
    print(f'The columns are {triad_tss_classification_continuous.columns}')
    plot_continuous(triad_tss_classification_continuous, threshold)