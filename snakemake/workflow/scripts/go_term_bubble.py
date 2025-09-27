import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# === Input files ===
broad_file = "../../../tmp/snakemake_intermediate_data/gProfiler_broad_results.csv"
narrow_file = "../../../tmp/snakemake_intermediate_data/gProfiler_narrow_results.csv"

# === Load data ===
broad = pd.read_csv(broad_file)
narrow = pd.read_csv(narrow_file)

# Add log-transformed values
broad["log_adj_p"] = -np.log10(broad["adjusted_p_value"])
narrow["log_adj_p"] = -np.log10(narrow["adjusted_p_value"])


# Filter for Biological Process (BP)
broad = broad[broad["source"] == "GO:BP"].sort_values(by=['log_adj_p'], ascending=False)
narrow = narrow[narrow["source"] == "GO:BP"].sort_values(by=['log_adj_p'], ascending=False)

# Keep top N terms by significance
topN = 5
broad_top = broad.head(5)
narrow_top = narrow.head(5)
print('broad')
print(broad_top)
print('narrow')
print(narrow_top)


# === Function to plot one dataset ===
def plot_bubble(df, title, outfile):
    plt.figure(figsize=(4, 6))
    
    # fixed x position
    x_vals = np.ones(len(df)) * 0.5
    
    scatter = plt.scatter(
        x=x_vals,
        y=df["term_name"],
        s=df["intersection_size"] ,
        c=df["log_adj_p"],
        cmap="viridis",
        alpha=0.8,
    )
    
    plt.gca().invert_yaxis() 
    plt.xticks([])
    plt.xlabel("")
    plt.title(title)
    
    # Add colorbar
    #cbar = plt.colorbar(scatter)
    #cbar.set_label("-log10(adjusted p-value)")
    
    plt.tight_layout()
    sns.despine()
    plt.savefig(outfile, dpi=900)
    plt.close()

# === Make two separate plots ===
plot_bubble(broad_top, "Broad", "broad_top5.png")
plot_bubble(narrow_top, "Narrow", "narrow_top5.png")
