#!/usr/bin/env python3
"""
Check how narrow-only and broad-only promoter genes are distributed across
expression categories (tissue-specific, widespread, null).

Inputs:
  - genes_only_narrow.txt
  - genes_only_broad.txt
  - cadenza_tissue_specific_list_z*.csv
  - gene correspondence table for ID mapping
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

# ========= INPUTS =========
genes_narrow_file = "../../../tmp/snakemake_intermediate_data/genes_only_narrow.txt"
genes_broad_file  = "../../../tmp/snakemake_intermediate_data/genes_only_broad.txt"
expression_file   = "../../../tmp/snakemake_intermediate_data/cadenza_tissue_specific_list_z1.468448433935507.csv"
conversion_table  = "../../../input_data/iwgsc_refseq_all_correspondances.csv"
outdir = "../../../tmp/snakemake_intermediate_data/results_promoter_vs_expression"
# ==========================

os.makedirs(outdir, exist_ok=True)

# --- Read promoter category lists ---
genes_narrow = pd.read_csv(genes_narrow_file, header=None, names=["gene"])
genes_narrow["gene"] = genes_narrow["gene"].astype(str).str.replace(r"\.\d+$", "", regex=True)
genes_narrow["promoter_class"] = "narrow_only"

# Read broad-only
genes_broad = pd.read_csv(genes_broad_file, header=None, names=["gene"])
genes_broad["gene"] = genes_broad["gene"].astype(str).str.replace(r"\.\d+$", "", regex=True)
genes_broad["promoter_class"] = "broad_only"

promoter_df = pd.concat([genes_narrow, genes_broad], ignore_index=True)
print(f"Loaded {len(genes_narrow)} narrow-only and {len(genes_broad)} broad-only genes")

# --- Read expression classification ---
expr = pd.read_csv(expression_file)

# Collapse multiple rows per gene into one expression_class
def collapse_categories(x):
    if (x == "null").any():
        return "null"
    elif (x == "widespread").any():
        return "widespread"
    else:
        return "tissue_specific"

expr_summary = expr.groupby("gene")["tissue"].apply(collapse_categories).reset_index()
expr_summary = expr_summary.rename(columns={"tissue": "expression_class"})

# --- Map gene IDs to the same version ---
def map_gene_ids(gene_df, correspondence_df_file_path,
                 from_version='v1.1', to_version='v2.1', gene_column='gene'):
    """Map gene IDs between IWGSC versions using correspondence table."""
    corr = pd.read_csv(correspondence_df_file_path, sep=None, engine="python")
    df = gene_df.copy()
    df[gene_column] = df[gene_column].astype(str).str.strip()
    for col in (from_version, to_version):
        corr[col] = corr[col].astype(str).str.strip()
    # Merge
    mapped = df.merge(corr[[from_version, to_version]],
                      left_on=gene_column, right_on=from_version, how="left")
    mapped = mapped.dropna(subset=[to_version])
    mapped[gene_column] = mapped[to_version]
    mapped = mapped.drop(columns=[from_version, to_version])
    return mapped

# Expression file is v1.1, promoter lists are v2.1 → convert expression to v2.1
expr_mapped = map_gene_ids(expr_summary, conversion_table,
                           from_version="v1.1", to_version="v2.1", gene_column="gene")

# --- Merge promoter classes with expression categories ---
merged = promoter_df.merge(expr_mapped, on="gene", how="left")

# --- Summary table ---
summary = (
    merged.groupby(["promoter_class", "expression_class"])
          .size()
          .reset_index(name="count")
)
summary["percent"] = summary.groupby("promoter_class")["count"].transform(
    lambda x: 100 * x / x.sum()
)

print("\nDistribution of promoter classes across expression categories:")

# Relabel 
ids = {'broad_only':'Broad',
      'narrow_only':  'Narrow'}
summary['promoter_class'] = summary['promoter_class'].replace(ids, regex=True)
#

ids = {'tissue_specific':'Tissue Specific',
      'widespread':  'Constitutive'}
summary['Expression Class'] = summary['expression_class'].replace(ids, regex=True)


print(summary)

custom_palette = {
    "Tissue Specific": "#005d5d",  # red
    "Constitutive": "#1192e8"
}
order_classes = ["Narrow", "Broad"]


# --- Plot stacked bar chart ---
plt.figure(figsize=(5.7, 5.7))
sns.barplot(data=summary, x="promoter_class", y="percent", hue="Expression Class", order=order_classes, palette=custom_palette)
plt.ylabel("Percentage of genes (%)", fontsize=20)
plt.xlabel("", fontsize=20)
plt.title("")
plt.tight_layout()
sns.despine()
plt.savefig(f"promoter_width_vs_expression.png", dpi=900)

plt.close()

# --- Save outputs ---
merged.to_csv(f"{outdir}/gene_promoter_vs_expression.csv", index=False)
summary.to_csv(f"{outdir}/summary_promoter_vs_expression.csv", index=False)

