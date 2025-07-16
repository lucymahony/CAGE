import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

# ---- CONFIG ----
results_dir = Path("/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/")
barplot_file = results_dir / "go_barplots.pdf"
heatmap_file = results_dir / "go_heatmap.png"

# ---- LOAD AND COMBINE GO RESULTS ----
dfs = []
tissue_label_map = {
    "SIS_SLE_SRO_SSP_genes_unique_SIS": "Immature_Spike",
    "SIS_SLE_SRO_SSP_genes_unique_SLE": "Leaf",
    "SIS_SLE_SRO_SSP_genes_unique_SRO": "Root",
    "SIS_SLE_SRO_SSP_genes_unique_SSP": "Spike"
}


for csv_file in results_dir.glob("GO_SIS_SLE_SRO_SSP*"):
    df = pd.read_csv(csv_file)
    if df.empty:
        continue
    tissue = csv_file.stem.replace("GO_", "").replace("_BP", "")
    df["tissue"] = tissue_label_map.get(tissue, tissue) 
    df = df.nsmallest(15, "classicFisher")
    df["term_label"] = df["Term"] + " (" + df["GO.ID"] + ")"
    df["logp"] = -np.log10(df["classicFisher"].replace(0, 1e-300))
    dfs.append(df)

if not dfs:
    raise ValueError("No non-empty *_BP.csv files found.")

all_go = pd.concat(dfs)

# ---- FACETED BARPLOT ----
sns.set(style="whitegrid")
g = sns.FacetGrid(all_go, col="tissue", col_wrap=3, sharex=False, height=5)
g.map_dataframe(sns.barplot, y="term_label", x="logp", palette="viridis")
g.set_titles("{col_name}")
g.set_axis_labels("-log10(p-value)", "")
for ax in g.axes.flatten():
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=8)
    ax.tick_params(axis='y', labelsize=8)
plt.tight_layout()
plt.savefig(barplot_file)
print(f"Saved barplots to {barplot_file}")

# ---- HEATMAP ----
heatmap_data = all_go.pivot_table(index="term_label", columns="tissue", values="logp", fill_value=0)
plt.figure(figsize=(12, 10))
sns.heatmap(heatmap_data, cmap="Blues", linewidths=0.5, linecolor='gray')
plt.title("GO Enrichment (-log10 p-value)")
plt.tight_layout()
plt.savefig(heatmap_file, dpi=300)
print(f"Saved heatmap to {heatmap_file}")
