# This script looks at the number of CAGE peaks conserved between the different tissues

# Import libraries
import pandas as pd
import pickle
from tqdm import tqdm
import numpy as np


cage_peaks = pd.read_csv('/Users/mahony/Documents/1.core_wheat/data/nanocage_intermediate_files/post_align_processing_only_unique/para_clustered_tc_2023-07-20_unique_only.csv',
                             sep=',', index_col=0)

print(cage_peaks.head())
print(cage_peaks.columns)
print(cage_peaks.shape)

# How many times is a cage peak within 50 bp of another cage peak on a different tissue?
# def count_close_peaks(df):
#     count = 0
#     for idx, row in tqdm(df.iterrows()):
#         close_peaks = df[(df['seqnames'] == row['seqnames']) &
#                          (df['group_name'] != row['group_name']) &
#                          (df['dominant_ctss'].between(row['dominant_ctss'] - 50, row['dominant_ctss'] + 50))].shape[0]
#         count += close_peaks
#     return count // 2  # Dividing by 2 because each pair will be counted twice

def count_close_peaks(df):
    # Sort the dataframe by seqnames and then by dominant_ctss
    df_sorted = df.sort_values(by=['seqnames', 'dominant_ctss'])

    # Helper function to count close peaks in each group
    def count_within_group(group):
        group['close_peak'] = group['dominant_ctss'].diff().between(1, 50) & (group['group_name'] != group['group_name'].shift())
        return group

    # Group by seqnames and apply the helper function
    df_sorted = df_sorted.groupby('seqnames').apply(count_within_group)
    close_peaks_count = df_sorted['close_peak'].sum()
    total_peaks = len(df_sorted)
    not_close_peaks = total_peaks - close_peaks_count
    print(f'Total peaks: {total_peaks}', f'Close peaks: {close_peaks_count}', f'Not close peaks: {not_close_peaks}', sep='\n')

result = count_close_peaks(cage_peaks)
print(result)