import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logomaker
from scipy.stats import fisher_exact   # you already have scipy
import seaborn as sns 
# ----------------------------
# Config (edit as needed)
# ----------------------------
FIMO_NARROW = "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/tmp/snakemake_intermediate_data/fimo_narrow/fimo.tsv"
FIMO_BROAD  = "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/tmp/snakemake_intermediate_data/fimo_broad/fimo.tsv"
FASTA_NARROW = "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/tmp/snakemake_intermediate_data/promoters_narrow.fa"    
FASTA_BROAD  = "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/tmp/snakemake_intermediate_data/promoters_broad.fa"     
MEME_LIB = (
    "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/"
    "scratch/repos/arabidopsis_circadian_binary_classicalML/"
    "statistical_distributions/input_matrixes/core_moftis_meme.txt"
)
OUTDIR = "motif_enrichment_output"
os.makedirs(OUTDIR, exist_ok=True)

# ----------------------------
# Helpers
# ----------------------------
def count_fasta_headers(fasta_path: str) -> int:
    """Count promoters via '>' headers (true denominators)."""
    with open(fasta_path) as fh:
        return sum(1 for line in fh if line.startswith(">"))

def benjamini_hochberg(pvals):
    """Manual Benjamini–Hochberg FDR (returns adjusted p-values in original order)."""
    p = np.asarray(pvals, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    adj = np.empty(n, dtype=float)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = ranked[i] * n / rank
        val = min(val, prev)  # monotonic
        adj[i] = val
        prev = val
    out = np.empty(n, dtype=float)
    out[order] = adj
    return out.tolist()

def load_fimo(path: str) -> pd.DataFrame:
    """Load FIMO table."""
    return pd.read_csv(path, sep="\t", comment="#")

def meme_parse(meme_file: str) -> dict:
    """Parse MEME motif file -> {motif_id: probability DataFrame (A,C,G,T)}."""
    motifs = {}
    with open(meme_file) as f:
        lines = f.readlines()

    i = 0
    L = len(lines)
    while i < L:
        line = lines[i].strip()
        if line.startswith("MOTIF"):
            parts = line.split()
            # handle 'MOTIF <id> [name]' form; take second token as id
            motif_id = parts[1]
            # advance until letter-probability line
            i += 1
            while i < L and not lines[i].startswith("letter-probability matrix"):
                i += 1
            if i >= L:
                break
            header = lines[i].strip()
            m = re.search(r"w=\s*(\d+)", header)
            if not m:
                raise ValueError(f"Could not find width (w=) in: {header}")
            width = int(m.group(1))
            i += 1
            rows = []
            for _ in range(width):
                vals = list(map(float, lines[i].split()))
                rows.append(vals)
                i += 1
            df = pd.DataFrame(rows, columns=["A", "C", "G", "T"])
            motifs[motif_id] = df
        else:
            i += 1
    return motifs

def plot_logo(prob_df: pd.DataFrame, motif_id: str, outdir: str, dpi=300):
    info_df = logomaker.transform_matrix(prob_df, from_type="probability", to_type="information")
    plt.figure(figsize=(max(3, info_df.shape[0]/2.5), 2.6))
    logomaker.Logo(info_df, color_scheme="classic")
    plt.title(motif_id)
    plt.ylabel("Bits")
    plt.xlabel("Position")
    plt.tight_layout()
    outpath = os.path.join(outdir, f"{motif_id}_logo.png")
    plt.savefig(outpath, dpi=dpi)
    plt.close()
    print(f"[logo] {motif_id} -> {outpath}")

# ----------------------------
# Main
# ----------------------------
def main():
    # Load FIMO
    fimo_n = load_fimo(FIMO_NARROW)
    fimo_b = load_fimo(FIMO_BROAD)

    # True denominators from FASTA headers
    n_prom = count_fasta_headers(FASTA_NARROW)   # e.g., 5504
    b_prom = count_fasta_headers(FASTA_BROAD)    # e.g., 101969
    print(f"Narrow promoters: {n_prom}")
    print(f"Broad promoters : {b_prom}")

    # Total hits per motif (for per-promoter averages)
    hits_n = fimo_n["motif_id"].value_counts()
    hits_b = fimo_b["motif_id"].value_counts()

    # Presence per promoter (≥1 hit) for Fisher’s exact (more interpretable)
    # sequence_name is the promoter identifier in FIMO
    pres_n = fimo_n.groupby("motif_id")["sequence_name"].nunique()
    pres_b = fimo_b.groupby("motif_id")["sequence_name"].nunique()

    motifs = sorted(set(hits_n.index).union(set(hits_b.index)))

    rows = []
    for m_id in motifs:
        # hits (for averages)
        h_n = int(hits_n.get(m_id, 0))
        h_b = int(hits_b.get(m_id, 0))

        # presence (for Fisher table)
        p_n = int(pres_n.get(m_id, 0))  # promoters (narrow) with ≥1 hit
        p_b = int(pres_b.get(m_id, 0))  # promoters (broad)  with ≥1 hit

        # contingency table on promoter presence:
        #         has motif   no motif
        # narrow     p_n      n_prom - p_n
        # broad      p_b      b_prom - p_b
        table = [[p_n, n_prom - p_n],
                 [p_b, b_prom - p_b]]
        odds, pval = fisher_exact(table, alternative="two-sided")

        avg_narrow = h_n / n_prom
        avg_broad  = h_b / b_prom

        rows.append({
            "motif_id": m_id,
            "narrow_hits": h_n,
            "broad_hits":  h_b,
            "narrow_promoters_with_motif": p_n,
            "broad_promoters_with_motif":  p_b,
            "avg_hits_per_promoter_narrow": avg_narrow,
            "avg_hits_per_promoter_broad":  avg_broad,
            "odds_ratio_presence": odds,
            "p_value": pval
        })

    df = pd.DataFrame(rows)
    rename_map = {
        "POL007.1": "BREd",
        "POL006.1":"BREu"}
    df["motif_id"] = df["motif_id"].replace(rename_map)

    # Manual BH-FDR
    df["FDR"] = benjamini_hochberg(df["p_value"].tolist())
 
    # Add promoter-level proportions
    df["prop_narrow"] = df["narrow_promoters_with_motif"] / n_prom
    df["prop_broad"]  = df["broad_promoters_with_motif"] / b_prom
 
    # Save full table
    out_tsv = os.path.join(OUTDIR, "motif_enrichment_results.tsv")
    df.sort_values("FDR").to_csv(out_tsv, sep="\t", index=False)
    print(f"[table] {out_tsv}")

    # Plot percentages for significant motifs
    sig = df[df["FDR"] < 0.05].copy()
 
    xlab = sig["motif_id"].tolist()
    x = np.arange(len(xlab))
    w = 0.4

    plt.figure(figsize=(6, 6))
    plt.bar(x - w/2, sig["prop_narrow"] * 100, width=w,
            label="Narrow", color="#ff7eb6")
    plt.bar(x + w/2, sig["prop_broad"] * 100, width=w,
            label="Broad", color="#8a3ffc")

    plt.xticks(x, xlab, ha="right")
    plt.ylabel("Promoters with ≥1 motif (%)", fontsize=14)
    plt.legend()
    plt.tight_layout()
    sns.despine()
    plot_path = os.path.join(OUTDIR, "percentage_promoters_sig.png")
    plt.savefig(plot_path, dpi=900)
    plt.close()
    plt.savefig('proportion.png')
if __name__ == "__main__":
    main()
