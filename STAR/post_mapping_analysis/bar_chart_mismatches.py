# This script takes the output from the script count_mismatches.sh which is the txt files such as IS1_Aligned.out_XM_histogram.txt
# and LE3_uniquelymapped_XM_histogram.txt and makes 8 bar chart, two for each tissue, one for UM and one for UO,
# showing the number of mismatches of the mapped reads


# Import libraries
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import sem

# Load in the data

def get_individual_hitogram_data(tissue, name):
    df = pd.read_csv(f"{tissue}_{name}_XM_histogram.txt", sep="\s+", header=None)
    df.columns = ['number_of_reads', 'number_of_mismatches']
    df['name'] = name
    # Add a column for tissue
    df['tissue'] = tissue
    df['tissue'] = df['tissue'].str[:2]
    df['number_of_reads'] = df['number_of_reads'] / 1000000 # Number of reads in millions
    return df


def get_all_data(name):
    # Names either ['Aligned.out', 'uniquelymapped']
    tissues = ['IS1', 'IS2', 'IS3', 'SP1', 'SP2', 'SP3', 'LE1', 'LE2', 'LE3', 'RO1', 'RO2', 'RO3']
    # Concatenate all the dataframes together with list comprehension
    df = pd.concat([get_individual_hitogram_data(tissue, name) for tissue in tissues])
    return df

def plot_results():
    title_map = {'IS': 'Immature Spike', 'SP': 'Spike', 'RO': 'Root', 'LE': 'Leaf'}

    # Create a 4x2 grid of bar charts
    fig, axes = plt.subplots(4, 2, figsize=(20, 10))
    axes = axes.flatten()

    for name in ['Aligned.out', 'uniquelymapped']:
        df = get_all_data(name)
        # Get the mean and standard error of the mean for each number of mismatches
        mean_df = df.groupby(['tissue', 'number_of_mismatches']).mean().reset_index()
        sem_df = df.groupby(['tissue', 'number_of_mismatches']).sem().reset_index()
        for i, tissue_type in enumerate(['IS', 'SP', 'RO', 'LE']):
            ax = axes[i] if name == 'Aligned.out' else axes[i + 4]
            sns.barplot(x='number_of_mismatches', y='number_of_reads', data = mean_df[mean_df['tissue'] == tissue_type], ax=ax, ci=None)
            # Add error bars for standard error
            x_vals = range(len(mean_df[mean_df['tissue'].str.startswith(tissue_type)]['number_of_mismatches']))
            y_vals = mean_df[mean_df['tissue'].str.startswith(tissue_type)]['number_of_reads']
            error_vals = sem_df[sem_df['tissue'].str.startswith(tissue_type)]['number_of_reads']
            ax.errorbar(x_vals, y_vals, yerr=error_vals, fmt='none', capsize=5, color='black')
            # Overlay the individual data points as gray dots
            sns.stripplot(x='number_of_mismatches', y='number_of_reads', data=df[df['tissue'].str.startswith(tissue_type)], ax=ax, color='gray', jitter=0.2, size=5)
            ax.set_title(title_map.get(tissue_type, tissue_type))
            ax.set_xlabel('Number of Mismatches')
            ax.set_ylabel('Million Reads')

    plt.tight_layout()
    plt.savefig('number_mismapped_reads.pdf')

plot_results()
