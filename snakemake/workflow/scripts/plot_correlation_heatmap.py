# Plots the output of the correlation_matrix funciton in cage_r_analysis.R as a heatmap 

import sys 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def plot_correlation_heatmap(input_file_path, output_file_path):
    """
    input_file_path (str): Path to the CSV file containing the correlation matrix.
    output_file_path (str): Path to save the heatmap image.
    """
    correlation_matrix = pd.read_csv(input_file_path, index_col=0)
    if not correlation_matrix.columns.equals(correlation_matrix.index):
            raise ValueError("The columns and index of the correlation matrix must match.")

    plt.figure(figsize=(10, 8))
    
    sns.heatmap(correlation_matrix, annot=False, fmt=".2f", cmap=sns.color_palette("ch:start=.2,rot=-.3", as_cmap=True), cbar=True,
                square=True, linewidths=0.5)
    plt.title("Correlation Heatmap")
    plt.tight_layout()
    plt.savefig(output_file_path, dpi =900)



if __name__ == "__main__":
    input_file_path = str(sys.argv[1])
    output_file_path = str(sys.argv[2])
    
    plot_correlation_heatmap(input_file_path, output_file_path)
    print('Script finished')

    # source ~/.bashrc 
    # mamba activate /hpc-home/mahony/miniforge3
    # Example command = python plot_correlation_heatmap.py ../../../intermediate_data/snakemake_intermediate_data/correlation_matrix.csv ../../../intermediate_data/snakemake_intermediate_data/correlation_matrix_seaborn_heatmap.png
