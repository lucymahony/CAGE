# This script takes number_reads_mapping.csv which is the output from calculate_number_times_reads_map.sh
# and generates 4 bar charts, one for each of the 4 tissues, showing the number of reads that map x number of times.

# Import libraries
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import sem

# Load in the data
number_reads_mapping = pd.read_csv("number_reads_mapping.csv", sep=",", header=0)
# Make a new column for tissue type which is the first two characters of the column tissue
number_reads_mapping['tissue_type'] = number_reads_mapping['tissue'].str[:2]
# Number of million reads

number_reads_mapping['number_of_reads'] = number_reads_mapping['number_of_reads'] / 1000000
# Inspect the data
print(number_reads_mapping.head())

# Calculate the mean value for each number_of_hits_in_reference within each tissue_type
mean_df = number_reads_mapping.groupby(['tissue_type', 'number_of_hits_in_reference'])['number_of_reads'].mean().reset_index()
sem_df = number_reads_mapping.groupby(['tissue_type', 'number_of_hits_in_reference'])['number_of_reads'].apply(sem).reset_index()

# Find the unique tissue_types
tissue_types = number_reads_mapping['tissue_type'].unique()
# Map tissue_type to descriptive titles
title_map = {'IS': 'Immature Spike', 'SP': 'Spike', 'RO': 'Root', 'LE': 'Leaf'}

# Create a 2x2 grid of bar charts
fig, axes = plt.subplots(2, 2, figsize=(10, 10))
axes = axes.flatten()

for i, tissue_type in enumerate(tissue_types):
    ax = axes[i]

    # Plot the bar chart using the mean values
    sns.barplot(x='number_of_hits_in_reference', y='number_of_reads',
                data=mean_df[mean_df['tissue_type'] == tissue_type], ax=ax, ci='se')
    # Add error bars for standard error
    x_vals = range(len(mean_df[mean_df['tissue_type'] == tissue_type]['number_of_hits_in_reference']))
    y_vals = mean_df[mean_df['tissue_type'] == tissue_type]['number_of_reads']
    error_vals = sem_df[sem_df['tissue_type'] == tissue_type]['number_of_reads']
    ax.errorbar(x_vals, y_vals, yerr=error_vals, fmt='none', capsize=5, color='black')

    # Overlay the individual data points as gray dots
    sns.stripplot(x='number_of_hits_in_reference', y='number_of_reads', data=number_reads_mapping[number_reads_mapping['tissue_type'] == tissue_type],
                  ax=ax, color='gray', jitter=0.2, size=5)

    ax.set_title(title_map.get(tissue_type, tissue_type))
    ax.set_xlabel('Number of Hits in Reference')
    ax.set_ylabel('Million Reads')

plt.tight_layout()
plt.savefig('number_reads_mapping.pdf')
