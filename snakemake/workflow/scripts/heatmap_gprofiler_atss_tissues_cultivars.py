
import pandas as pd
import numpy as np
import os

# ---- Config ----
input_dir = "../../../results/Alternative_TSS_genes_tissues_varieties"
files = {
    "Spike": "gProfiler_SP_unique.csv",
    "ImmatureSpike": "gProfiler_IS_unique.csv",
    "Root": "gProfiler_RO_unique.csv",
    "Leaf": "gProfiler_LE_unique.csv"
}
padj_cutoff = 0.05
top_terms_per_tissue=5

# ---- Collect all results ----
rows = []

for tissue, filename in files.items():
    print(tissue)
    filepath = os.path.join(input_dir, filename)
    df = pd.read_csv(filepath)
    print(df.head())

    # filter significant
    df = df[df["adjusted_p_value"] <= padj_cutoff].copy()
    # Filter BP
    #df = df[df['source'] == 'GO:BP']
    if df.empty:
        print(f"No significant terms in {tissue}")
        continue

    # add log10 transform
    df["-log10(padj)"] = -np.log10(df["adjusted_p_value"].replace(0, 1e-300))
    df = df.nsmallest(top_terms_per_tissue, "adjusted_p_value")

    for _, row in df.iterrows():
        rows.append({
            "Tissue": tissue,
            "Source": row["source"],
            "GO Term": row["term_id"],
            "Description": row["term_name"],
            "-log10(padj)": round(row["-log10(padj)"], 2,)
        })

# ---- Build table ----
results_table = pd.DataFrame(rows)

# Save as CSV
results_table.to_csv("GO_summary_table.csv", index=False)


print("Wrote GO_summary_table.csv and GO_summary_table.md")
print(results_table)
