# This script takes number_reads_mapping.csv which is the output from calculate_number_times_reads_map.sh
# and generates 4 bar charts, one for each of the 4 tissues, showing the number of reads that map x number of times.

# Import libraries
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load in the data
number_reads_mapping = pd.read_csv("number_reads_mapping.csv", sep=",", header=0)
# Make a new column for tissue type which is the first two characters of the column tissue
number_reads_mapping['tissue_type'] = number_reads_mapping['tissue'].str[:2]

print(number_reads_mapping.head())
