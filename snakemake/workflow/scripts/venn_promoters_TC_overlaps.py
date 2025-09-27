import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles
import venn

# ---- Load data ----
file = "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/tmp/snakemake_intermediate_data/tmp_overlap_matrix/ALL_with_overlap_matrix.csv"
output_file_path = "../../../results/Alternative_TSS_genes_tissues_varieties"
cols = [
    "id","seqnames","start","end","score","strand",
    "dominant_start","dominant_end","annotation","nearest_gene",
    "extra1","extra2",  # junk columns
    "CL","CR","SIS","SLE","SRO","SSP"
]

df = pd.read_csv(file, names=cols, header=0, low_memory=False)
df = df.drop(columns=["extra1","extra2"])
df = df[df["annotation"] == "promoter"]

# Force numeric conversion for counts
num_cols = ["start","end","score","dominant_start","dominant_end",
            "CL","CR","SIS","SLE","SRO","SSP"]
for col in num_cols:
    df[col] = pd.to_numeric(df[col], errors="coerce")

# ---- Mapping dictionaries ----
pretty_names = {
    "CL": "Cadenza leaf",
    "CR": "Cadenza root",
    "SIS": "Fielder immature spike",
    "SLE": "Fielder leaf",
    "SRO": "Fielder root",
    "SSP": "Fielder spike"
}

# Colors match matplotlibs tab§10
colors = {
    "Cadenza leaf": "#1f77b4",     # tab10 blue
    "Cadenza root": "#ff7f0e",     # tab10 orange
    "Fielder immature spike": "#9ecae1",  # pastel blue
    "Fielder leaf": "#fcae91",     # pastel red
    "Fielder root": "#f2b6f5",     # pastel purple
    "Fielder spike": "#9fe5eb"     # pastel teal
}

# ---- Convert to sets of genes ----
def get_genes(df, col):
    return set(df.loc[df[col] == 1, "nearest_gene"].dropna().astype(str)) - {""}

sets = {k: get_genes(df, k) for k in ["CL","CR","SIS","SLE","SRO","SSP"]}

# ---- Unique gene lists (per tissue) ----
for key, genes in sets.items():
    filename = f"{pretty_names[key].replace(' ', '_')}_unique_genes.txt"
    with open(filename, "w") as f:
        for g in sorted(genes):
            f.write(g + "\n")
    print(f"Wrote {len(genes)} genes to {filename}")

# ---- Unique to Cadenza vs Fielder ----
cadenza_union = sets["CL"].union(sets["CR"])
fielder_union = sets["SIS"].union(sets["SLE"], sets["SRO"], sets["SSP"])

cadenza_unique = cadenza_union - fielder_union
fielder_unique = fielder_union - cadenza_union

with open("Cadenza_unique_genes.txt", "w") as f:
    for g in sorted(cadenza_unique):
        f.write(g + "\n")

with open("Fielder_unique_genes.txt", "w") as f:
    for g in sorted(fielder_unique):
        f.write(g + "\n")

print(f"Wrote {len(cadenza_unique)} unique Cadenza genes")
print(f"Wrote {len(fielder_unique)} unique Fielder genes")

# ---- Venn plotting helpers ----
def venn_plot(sets_dict, filename, title):
    """Plot with venn (supports >2 sets)."""
    print(f"\n{title}")
    for name, s in sets_dict.items():
        print(f"  {pretty_names.get(name, name)}: {len(s)} genes")

    if all(len(s) == 0 for s in sets_dict.values()):
        print("  Skipping plot: all sets empty")
        return

    pretty_sets = {pretty_names[k]: v for k, v in sets_dict.items()}
    plt.figure(figsize=(6,6))
    venn.venn(pretty_sets, cmap="tab10")
    plt.title(title)
    plt.savefig(filename, dpi=900)
    plt.close()

def matplot_ven(set1, set2, labels, filename, title):
    """Plot 2-way venn with matplotlib_venn, styled with hex colours."""
    plt.figure(figsize=(6,6))
    v = venn2([set1, set2], set_labels=labels)

    # Colour the circles
    for i, label in enumerate(labels):
        patch_id = "10" if i == 0 else "01"
        if v.get_patch_by_id(patch_id):
            v.get_patch_by_id(patch_id).set_color(colors[label])
            v.get_patch_by_id(patch_id).set_alpha(0.6)

    # Add coloured outlines
    c = venn2_circles([set1, set2])
    for i, label in enumerate(labels):
        c[i].set_edgecolor(colors[label])
        c[i].set_linewidth(2)

    plt.title(title)
    plt.savefig(filename, dpi=900)
    plt.close()

# ---- Example plots ----
# Pairwise (matplotlib_venn)
matplot_ven(sets["CL"], sets["SLE"],
            labels=(pretty_names["CL"], pretty_names["SLE"]),
            filename=f"{output_file_path}/matplot_CL_vs_SLE.png",
            title="")

matplot_ven(sets["CR"], sets["SRO"],
            labels=(pretty_names["CR"], pretty_names["SRO"]),
            filename=f"{output_file_path}/matplot_CR_vs_SRO.png",
            title="")

matplot_ven(sets["CL"], sets["CR"],
            labels=(pretty_names["CL"], pretty_names["CR"]),
            filename=f"{output_file_path}/matplot_CL_vs_CR.png",
            title="")

matplot_ven(sets["SLE"], sets["SRO"],
            labels=(pretty_names["SLE"], pretty_names["SRO"]),
            filename=f"{output_file_path}/matplot_SLE_vs_SRO.png",
            title="")

# Multi-way (venn lib)

venn_plot({"SIS": sets["SIS"], "SLE": sets["SLE"], "SRO": sets["SRO"], "SSP": sets["SSP"]},
          f"{output_file_path}/venn_Fielder_tissues.png", "")
