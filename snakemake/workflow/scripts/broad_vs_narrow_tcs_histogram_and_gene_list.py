# Generates a histogram of the widths
# Note set  the width cut off
# Note decide to filter or not 


import matplotlib.pyplot as plt
import pandas as pd 
import seaborn as sns

csv_file = "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/tmp/snakemake_intermediate_data/no_merge_reps/finalised_custom_annotations/ALL_TC_merged_output.csv"
width_cutoff = 10   
df = pd.read_csv(csv_file, index_col=False)

df["score"] = pd.to_numeric(df["score"], errors="coerce")
df["start"] = pd.to_numeric(df["start"], errors="coerce")
df["end"] = pd.to_numeric(df["end"], errors="coerce")



df["width"] = df["end"] - df["start"] + 1
# Basic stats
min_width = df["width"].min()
median_width = df["width"].median()
mean_width = df['width'].mean()
max_width = df["width"].max()


print(f"Width statistics (bp):")
print(f"  Min: {min_width}")
print(f"  Median: {median_width}")
print(f"  Mean: {mean_width}")
print(f"  Max: {max_width}")

less_than = df[df["width"] < width_cutoff]["width"]
greater_equal = df[df["width"] >= width_cutoff]["width"]

# Print results
total = len(df)
print(f"Total tag clusters (score > 0.5): {total}")
print(f"Widths < {width_cutoff} bp: {len(less_than)} ({len(less_than)/total*100:.2f}%)")
print(f"Widths ≥ {width_cutoff} bp: {len(greater_equal)} ({len(greater_equal)/total*100:.2f}%)")

# Plot histogram
plt.hist([less_than, greater_equal],
         bins=100, range=(0,100),
         stacked=True,
         color=['#ff7eb6', '#be95ff'], # "#9f1853", "#6929c4"
         label=[f"< {width_cutoff} bp", f"≥ {width_cutoff} bp"],)

plt.xlabel("Tag Cluster Width (bp)", fontsize=20)
plt.ylabel("Frequency", fontsize=20)
plt.title(f"")
plt.tight_layout()
sns.despine()
plt.savefig("tag_cluster_width_distribution.png", dpi=900)



##### List geness in narrow and broad category
outdir = "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/tmp/snakemake_intermediate_data/"

df["nearest_gene"] = df["nearest_gene"].fillna("")

# Group by gene and classify
gene_groups = df[df["nearest_gene"] != ""].groupby("nearest_gene")["width"].apply(list)


genes_only_narrow = []
genes_only_broad = []
genes_mixed = []

for gene, widths in gene_groups.items():
    if gene == "":
        continue
    if all(w < width_cutoff for w in widths):
        genes_only_narrow.append(gene)
    elif all(w >= width_cutoff for w in widths):
        genes_only_broad.append(gene)
    else:
        genes_mixed.append(gene)

# Write outputs
with open(f"{outdir}/genes_only_narrow.txt", "w") as f:
    for g in genes_only_narrow:
        f.write(g + "\n")

with open(f"{outdir}/genes_only_broad.txt", "w") as f:
    for g in genes_only_broad:
        f.write(g + "\n")

with open(f"{outdir}/genes_mixed.txt", "w") as f:
    for g in genes_mixed:
        f.write(g + "\n")

# Print summary
print(f"Genes with only narrow clusters (<{width_cutoff} bp): {len(genes_only_narrow)}")
print(f"Genes with only broad clusters (≥{width_cutoff} bp): {len(genes_only_broad)}")
print(f"Genes with both narrow and broad clusters: {len(genes_mixed)}")