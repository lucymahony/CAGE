# Plots the IQR widths from the ../../../intermediate_data/snakemake_intermediate_data/shared_fielder_cadenza_all_at_least_6_overlaps_max_score.bed
# Chr1A	1245	1270	.	24.1729202084705	+	1249	1250
# Chr1A	5235	5265	.	55.3547183044327	+	5239	5240
# Chr1A	6289	6344	.	85.008305053579	+	6330	6331


import sys 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def concatenate_multiple_TC_width_files(input_file_directory, output_file_path):
    list_of_input_files = ['SSP3_TC_width.txt',
                        'SSP2_TC_width.txt',
                        'SSP1_TC_width.txt',
                        'SRO3_TC_width.txt',
                        'SRO2_TC_width.txt',
                        'SRO1_TC_width.txt',
                        'SLE3_TC_width.txt',
                        'SLE2_TC_width.txt',
                        'SLE1_TC_width.txt',
                        'SIS3_TC_width.txt',
                        'SIS2_TC_width.txt',
                        'SIS1_TC_width.txt',
                        'CR5_TC_width.txt',
                        'CR4_TC_width.txt',
                        'CR3_TC_width.txt',
                        'CR2_TC_width.txt',
                        'CR1_TC_width.txt',
                        'CL5_TC_width.txt',
                        'CL4_TC_width.txt',
                        'CL3_TC_width.txt',
                        'CL2_TC_width.txt',
                        'CL1_TC_width.txt']

    # Read all the files in as pd dataframes and the concatenate them
    all_files = []
    for file in list_of_input_files:
        file_path = input_file_directory + file
        csv = pd.read_csv(file_path, sep="	", skiprows=1)
        csv.columns = ['seqnames', 'start', 'end', 'width', 'strand', 'score', 'dominant_ctss.seqnames', 'dominant_ctss.pos', 'dominant_ctss.strand', 'dominant_ctss.score', 'nr_ctss']
        print(csv.head())
        all_files.append(csv) # Skip the header
        print('File read: ', file_path)
    concatenated_files = pd.concat(all_files)
    # Save concatenated files
    file_path = output_file_path + '.csv'
    concatenated_files.to_csv(file_path, sep=',', index=False)
    print('Concatenated files saved to: ', file_path)
    return concatenated_files


def plot_IQR_widths_histogram(tc_widths_concatenated, output_file_path):
    """
    input_file_path (df): Dataframe with the IQR widths of the tag clusters
    output_file_path (str): Path to save the distribution plot 
    """

    print('Mean width: ', tc_widths_concatenated['width'].mean())
    print('Std width: ', tc_widths_concatenated['width'].std())

    # Remove rows with width > 150 for plotting 
    tc_widths_concatenated = tc_widths_concatenated[tc_widths_concatenated['width'] < 150]
    # Plot the width histogram
    sns.histplot(tc_widths_concatenated['width'], bins=150, kde=True)
    plt.xlabel('Width')
    plt.ylabel('Frequency')
    plt.title('Distribution of IQR widths')
    file_name = output_file_path +  '_histplot.pdf'
    plt.savefig(file_name, dpi = 900)
    print('Plot saved to: ', file_name)

def plot_IQR_widths_kde(tc_widths_concatenated, output_file_path):
    # Remove rows with width > 150 for plotting 
    tc_widths_concatenated = tc_widths_concatenated[tc_widths_concatenated['width'] < 150]
    # Plot the width with kde plot 
    plt.figure(figsize=(5.7, 5.7))
    sns.kdeplot(data=tc_widths_concatenated, x='width', fill=True, color='#3ddbd9', alpha=1, linewidth=0)
    plt.xlabel('Tag Cluster Width (bp)', fontsize=11)
    plt.ylabel('Frequency', fontsize=11)
    plt.title('')
    plt.xlim(0, 150)
    file_name = output_file_path +  '_kdeplot.pdf'
    sns.despine()
    plt.savefig(file_name, dpi = 900)
    print('Plot saved to: ', file_name)
    



if __name__ == "__main__":
    input_file_directory = str(sys.argv[1]) 
    output_file_path = str(sys.argv[2])

    #concatenated_files = concatenate_multiple_TC_width_files(input_file_directory, output_file_path)

    concatenated_files = pd.read_csv(output_file_path + '.csv', sep=',', skiprows=1)
    concatenated_files.columns = ['seqnames', 'start', 'end', 'width', 'strand', 'score', 'dominant_ctss.seqnames', 'dominant_ctss.pos', 'dominant_ctss.strand', 'dominant_ctss.score', 'nr_ctss']


    print('Concatenated files: ', concatenated_files.head())
    plot_IQR_widths_kde(concatenated_files, output_file_path)
    
    print('Script finished')

    # source ~/.bashrc 
    # mamba activate /hpc-home/mahony/miniforge3
    # Example command = python plot_iqr_widths_distribution.py ../../../intermediate_data/snakemake_intermediate_data/ ../../../intermediate_data/snakemake_intermediate_data/shared_fielder_cadenza_all_TC_widths
