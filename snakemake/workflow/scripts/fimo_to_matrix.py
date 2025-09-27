import pandas as pd
import re

# === INPUT FILES ===
fimo_path = "fimo_out/fimo.tsv"
fasta_path = "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/tmp/snakemake_intermediate_data/no_merge_reps/finalised_custom_annotations/high_conf_genes_promoters.fa"


# === Load FIMO hits ===
fimo = pd.read_csv(fimo_path, sep="\t")
header_map = {}  # (chr, start, end, strand) → gene
with open(fasta_path) as f:
    for line in f:
        if line.startswith(">"):
            line = line.strip()[1:]  # remove ">"
            # e.g., TraesCS1A03G0004200::Chr1A:1161066-1161216(-)
            m = re.match(r"(?P<gene>[^:]+)::(?P<chr>[^:]+):(?P<start>\d+)-(?P<end>\d+)\((?P<strand>[+-])\)", line)
            if m:
                key = (
                    m["chr"],
                    int(m["start"]),
                    int(m["end"]),
                    m["strand"]
                )
                header_map[key] = m["gene"]

# === Step 3: Map FIMO hits to gene names ===
def find_gene(chr, start, end, strand):
    for (g_chr, g_start, g_end, g_strand), gene in header_map.items():
        if g_chr == chr and g_strand == strand and g_start <= start <= g_end:
            return gene
    return None

fimo["gene"] = fimo.apply(
    lambda row: find_gene(row["sequence_name"], row["start"], row["stop"], row["strand"]),
    axis=1
)

# === Step 4: Build gene × motif count matrix ===
motif_matrix = (
    fimo.dropna(subset=["gene"])
        .groupby(["gene", "motif_id"])
        .size()
        .unstack(fill_value=0)
        .reset_index()
)

# === Step 5: Save output ===
motif_matrix.to_csv("motif_count_matrix.csv", index=False)
print("✔ Saved motif count matrix to: motif_count_matrix.csv")