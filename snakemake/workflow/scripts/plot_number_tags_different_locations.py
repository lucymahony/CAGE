
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd 

# Data 1 
samples = [
    "L_cadenza", "R_cadenza", "RO_fielder", "LE_fielder", "SP_fielder", 
    "IS_fielder"
]
high_confidence = [27768, 30301, 27339, 16268, 26524, 26676]
low_confidence = [2474, 2979, 2273, 1138, 2443, 2958]
enhancers = [142735, 177026, 139620, 58298, 139021, 177157]
cis_nats = [50605, 55886, 46658, 24851, 46688, 46464]

df = pd.DataFrame({
    'high_confidence': high_confidence,
    'low_confidence': low_confidence,
    'enhancers': enhancers,
    'cis_nats': cis_nats
}, index=samples)



# Bar chart with high_confidence, low_confidence, enhancers and cis_nats as the x axis, with error bars as the standard deviation between the samples
# Y values are the mean of the samples
x = np.array(['High\nConfidence\nGenes', 'Low\nConfidence\nGenes', 'Putative\nEnhancers', 'Cis-NATs'])
y = np.array([df['high_confidence'].mean(), df['low_confidence'].mean(), df['enhancers'].mean(), df['cis_nats'].mean()])
yerr = np.array([df['high_confidence'].std(), df['low_confidence'].std(), df['enhancers'].std(), df['cis_nats'].std()])

plt.figure(figsize=(5.7, 5.7))
colours = ['#ff8389', '#ff7eb6', '#be95ff', '#78a9ff']
plt.bar(x, y, yerr=yerr, capsize=5, color=colours)
plt.ylabel("Number of Tag Clusters")
plt.gca().set_yscale('log')
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.savefig("number_tags_different_locations.png", dpi=900)


# Dataset 2
samples = [
    "Overlap in every tissue in fielder and cadenza", 
    "Only 1 tissue", "Total number"
]
high_confidence = [6732, 11967, 48270]
low_confidence = [455, 1992, 5381]
enhancers = [46655, 11906, 472576]
cis_nats = [20046, 7274, 64685]
# Order is Shared, Unique, total 

# Total, Shared, Unique for high_confidence, low_confidence, enhancers and cis_nats
df = pd.DataFrame({
    'high_confidence': high_confidence,
    'low_confidence': low_confidence,
    'enhancers': enhancers,
    'cis_nats': cis_nats
}, index=samples)

# Bar chart of shared, unique and total for high_confidence, low_confidence, enhancers and cis_nats
x = np.array(['Shared\nHigh\nConfidence\nGenes', 'Unique\nHigh\nConfidence\nGenes', 'Total\nHigh\nConfidence\nGenes', 'Shared\nLow\nConfidence\nGenes', 'Unique\nLow\nConfidence\nGenes', 'Total\nLow\nConfidence\nGenes', 'Shared\nPutative\nEnhancers', 'Unique\nPutative\nEnhancers', 'Total\nPutative\nEnhancers', 'Shared\nCis-NATs', 'Unique\nCis-NATs', 'Total\nCis-NATs'])
y = np.array([
    df['high_confidence'][0], df['high_confidence'][1], df['high_confidence'][2],
    df['low_confidence'][0], df['low_confidence'][1], df['low_confidence'][2],
    df['enhancers'][0], df['enhancers'][1], df['enhancers'][2],
    df['cis_nats'][0], df['cis_nats'][1], df['cis_nats'][2]
])

plt.figure(figsize=(5.7, 5.7))
plt.xticks(fontsize=3)
colours = ['#ff8389', '#fa4d56', '#da1e28', '#ff7eb6', '#ee5396', '#d02670', '#be95ff', '#a56eff', '#8a3ffc', '#78a9ff', '#4589ff', '#0f62fe']
plt.bar(x, y, color=colours)
plt.ylabel("Number of Tag Clusters")
plt.gca().set_yscale('log')
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.savefig("number_tags_shared_unique_high_confidence.png", dpi=900)

