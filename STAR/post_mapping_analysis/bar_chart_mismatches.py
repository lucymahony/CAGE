# This script takes the output from the script count_mismatches.sh which is the txt files such as IS1_Aligned.out_XM_histogram.txt
# and LE3_uniquelymapped_XM_histogram.txt and makes 8 bar chart, two for each tissue, one for UM and one for UO,
# showing the number of mismatches of the mapped reads


# Import libraries
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import sem

# Load in the data

def get_data(tissue, name):
    df = pd.read_csv(f"{tissue}_{name}_XM_histogram.txt", sep=" ", header=None)
    df.columns = ['number_of_reads', 'number_of_mismatches']
    df['tissue'] = tissue
    df['name'] = name
    return df

IS1_UM = get_data('IS1', 'Aligned.out')

print(IS1_UM.head())