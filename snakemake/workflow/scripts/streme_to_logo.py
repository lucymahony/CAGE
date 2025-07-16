import os
import sys
import pandas as pd
import logomaker
import matplotlib.pyplot as plt

def parse_meme_file(meme_file):
    motifs = {}
    with open(meme_file) as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("MOTIF"):
            motif_id = line.strip().split()[1]
            # Find letter-probability matrix
            while not lines[i].startswith("letter-probability matrix"):
                i += 1
            header = lines[i]
            tokens = header.strip().split()
            length = int(tokens[tokens.index("w=") + 1])
            i += 1
            matrix = []
            for _ in range(length):
                values = list(map(float, lines[i].strip().split()))
                matrix.append(values)
                i += 1
            df = pd.DataFrame(matrix, columns=["A", "C", "G", "T"])
            motifs[motif_id] = df
        else:
            i += 1
    return motifs

def plot_logo_info_content(prob_df, motif_id, output_dir, dpi=600):
    # Convert probabilities to information content
    info_df = logomaker.transform_matrix(prob_df, from_type='probability', to_type='information')

    # Plot
    plt.figure(figsize=(info_df.shape[0] / 2.5, 2.5))
    logo = logomaker.Logo(info_df, font_name='Arial', color_scheme='classic')

    logo.ax.set_title(motif_id)
    logo.ax.set_ylabel("Bits")
    logo.ax.set_xlabel("Position")
    logo.ax.set_ylim([0, 2])  # Optional: cap at 2 bits
    logo.style_spines(visible=False)
    logo.style_spines(spines=['left', 'bottom'], visible=True)
    plt.tight_layout()
    output_path = os.path.join(output_dir, f"{motif_id}_info_logo.png")
    plt.savefig(output_path, dpi=dpi)
    plt.close()
    print(f"Saved: {output_path}")

def main(meme_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    motifs = parse_meme_file(meme_file)
    for motif_id, prob_df in motifs.items():
        plot_logo_info_content(prob_df, motif_id, output_dir)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python streme_to_logo.py <streme.txt> <output_dir>")
        sys.exit(1)
    
    meme_file = sys.argv[1]
    output_dir = sys.argv[2]
    main(meme_file, output_dir)
