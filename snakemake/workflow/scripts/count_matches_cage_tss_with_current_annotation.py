import pandas as pd
import argparse
import os


def parse_arguments():
    parser = argparse.ArgumentParser(description="Count matched CAGE TSSs against wheat gene annotations.")
    parser.add_argument('--cage_csv', default='/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/tmp/snakemake_intermediate_data/no_merge_reps/finalised_custom_annotations/ALL_TC_merged_output.csv')
    parser.add_argument('--gff3_hc', default='/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/input_data/chinese_spring_genome_data/iwgsc_refseqv2.1_annotation_200916_HC.gff3')
    parser.add_argument('--gff3_lc', default='/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/input_data/chinese_spring_genome_data/iwgsc_refseqv2.1_annotation_200916_LC.gff3')
    parser.add_argument('--tss_window', type=int, default=100, help="Window size (bp) around annotated TSS for matching")
    parser.add_argument('--fasta_file', default='/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Wheat_expression_prediction/intermediate_data/promoter_seqs_150_0.fa')
    parser.add_argument('--output_file', default='/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/tmp/snakemake_intermediate_data/no_merge_reps/finalised_custom_annotations/high_conf_genes_promoters.fa')
    parser.add_argument('--summary_results_file', default='/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/tmp/snakemake_intermediate_data/no_merge_reps/finalised_custom_annotations/ALL_TC_merged_SUMMARY_STATS.txt')
    return parser.parse_args()


def load_cage_tss(csv_file):
    colnames = [
        "seqnames", "start", "end", "score", "strand",
        "dominant_start", "dominant_end", "annotation",
        "nearest_gene", "nearest_gene_2"
    ]

    # Read malformed rows safely with explicit names
    df = pd.read_csv(
        csv_file,
        header=0,
        names=colnames,
        quotechar='"',
        skipinitialspace=True,
        engine='python'
    )

    # Drop rows with missing gene info
    df = df[df['nearest_gene'].notnull() & (df['nearest_gene'] != '')]
    df['nearest_gene'] = df['nearest_gene'].str.replace(r'\.\d+$', '', regex=True)


    # Convert score and dominant_start
    df['score'] = df['score'].astype(float)
    #df = df[df['score'] > 0.5]

    df['tss'] = df['dominant_start'].astype(int)

    df = df[['seqnames', 'strand', 'tss', 'nearest_gene']]
    df = df.rename(columns={
        'seqnames': 'chrom',
        'nearest_gene': 'gene_id'
    })

    return df


def parse_gff3_tss(gff3_file):
    rows = []
    with open(gff3_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            if len(parts) != 9:
                continue
            chrom, source, feature, start, end, score, strand, phase, attributes = parts
            if feature not in ('mRNA', 'transcript'):
                continue
            start, end = int(start), int(end)
            gene_id = None
            for attr in attributes.split(';'):
                if attr.startswith('Parent='):
                    gene_id = attr.split('=')[1]
                    break
                elif attr.startswith('gene_id='):
                    gene_id = attr.split('=')[1]
            if gene_id:
                tss = start if strand == '+' else end
                rows.append({'chrom': chrom, 'tss': tss, 'strand': strand, 'gene_id': gene_id})
    return pd.DataFrame(rows)


def count_matches(cage_df, annot_df, window, save_matched_genes=False):
    # Merge by gene ID and chromosome and strand
    merged = pd.merge(cage_df, annot_df, on=['gene_id', 'chrom', 'strand'], suffixes=('_cage', '_annot'))

    # Compute absolute TSS distance
    merged['tss_diff'] = (merged['tss_cage'] - merged['tss_annot']).abs()

    # Mark matches within window
    matched = merged[merged['tss_diff'] <= window]
    unmatched = merged[merged['tss_diff'] > window]

    total_genes_with_cage = cage_df['gene_id'].nunique()
    matched_genes = matched['gene_id'].nunique()

    if save_matched_genes:
        # matched
        matched_gene_ids = matched['gene_id'].unique()
        output_file=f'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/results/genes_hc_tss_validated_by_cage_within_{window}bp.txt'
        with open(output_file, 'w') as f:
            for gene_id in matched_gene_ids:
                f.write(f"{gene_id}\n")
        # unmatched
        unmatched_gene_ids = unmatched['gene_id'].unique()
        output_file=f'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/results/genes_hc_tss_not_validated_by_cage_within_{window}bp.txt'
        with open(output_file, 'w') as f:
            for gene_id in unmatched_gene_ids:
                f.write(f"{gene_id}\n")

    return total_genes_with_cage, matched_genes

def extract_promoters_from_gene_list(gene_list_path, fasta_path, output_path):
    # Load list of gene IDs
    with open(gene_list_path, 'r') as f:
        gene_ids = set(line.strip() for line in f if line.strip())

    # Prepare output
    with open(fasta_path, 'r') as fasta_in, open(output_path, 'w') as fasta_out:
        write_record = False
        header = ''
        seq = ''

        for line in fasta_in:
            line = line.rstrip()
            if line.startswith(">"):
                # Save the previous record if it matched
                if write_record and header and seq:
                    fasta_out.write(f"{header}\n{seq}\n")
                # Reset
                seq = ''
                header = line
                # Get gene ID from FASTA header (before ::)
                gene_id = line[1:].split("::")[0]
                write_record = gene_id in gene_ids
            else:
                if write_record:
                    seq += line

        # Don't forget the last record
        if write_record and header and seq:
            fasta_out.write(f"{header}\n{seq}\n")

    print(f"Promoter sequences for {len(gene_ids)} genes written to: {output_path}")

def main():
    args = parse_arguments()

    with open(args.summary_results_file, 'w') as summary:
        def log(msg):
            print(msg)
            summary.write(str(msg) + "\n")

        log(f"Loading CAGE data from:\n  {args.cage_csv}")
        cage_df = load_cage_tss(args.cage_csv)

        log(f"Parsing annotation from:\n  {args.gff3_hc}")
        annot_df = parse_gff3_tss(args.gff3_hc)
        log(str(annot_df.head()))

        log("\nResults:")

        window = 1
        total, matched = count_matches(cage_df, annot_df, window)
        log(f"Total genes with CAGE TSSs: {total}")
        log(f"Genes with matching annotated TSSs (±{window} bp): {matched} = {matched/total*100:.2f}%")

        window = 10
        total, matched = count_matches(cage_df, annot_df, window, save_matched_genes=True)
        log(f"Genes with matching annotated TSSs (±{window} bp): {matched} = {matched/total*100:.2f}%")

        validated_gene_list = f'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/results/genes_hc_tss_validated_by_cage_within_{window}bp.txt'


        extract_promoters_from_gene_list(
            validated_gene_list,
            args.fasta_file,
            args.output_file
        )

        log(f"\nSummary results written to: {args.summary_results_file}")

if __name__ == "__main__":
    main()
