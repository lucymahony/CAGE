#!/bin/bash 
#SBATCH -p ei-medium
#SBATCH -o compare.out
#SBATCH -c 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk
#SBATCH --mem=12G 
#SBATCH --time=0-00:30:00

source package 4028d6e4-21a8-45ec-8545-90e4ed7e1a64  # BedTools2 2.30.0

# Input files
cage_tag_clusters=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/no_merge_reps/TC_GR_CL5_no_merge_distclu_md50_single_5_myCAGEset_my_annot.csv
reference_annotation=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/input_data/chinese_spring_genome_data/iwgsc_refseqv2.1_annotation_200916_HC.gff3

# Output files
directory=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/
cage_tag_clusters_bed="${cage_tag_clusters%.csv}_tc_regions.bed"
cage_tag_dominant_ctss_bed="${cage_tag_clusters%.csv}_dominant_ctss.bed"
reference_annotation_tss_bed="${reference_annotation%.gff3}_tss.bed"
reference_annotation_tss_10_bed="${reference_annotation%.gff3}_tss_10bp_window.bed"

# Sorted versions
sorted_dominant_ctss="${directory}sorted_dominant_ctss.bed"
sorted_tc_regions="${directory}sorted_tc_regions.bed"
sorted_annotation_tss="${directory}sorted_annotation_tss.bed"
sorted_annotation_tss_10="${directory}sorted_annotation_tss_10bp_window.bed"

echo "Creating BED files..."

# Extract dominant CTSS
awk -F',' 'NR>1 {
    gsub(/"/, "", $7); gsub(/"/, "", $9);
    OFS="\t"; print $7, $8-1, $8, "TC", ".", $9
}' "$cage_tag_clusters" > "$cage_tag_dominant_ctss_bed"

# Extract TC regions
awk -F',' 'NR>1 {
    gsub(/"/, "", $1); gsub(/"/, "", $2); gsub(/"/, "", $3); gsub(/"/, "", $5);
    OFS="\t"; print $1, $2, $3, "TC_region", ".", $5
}' "$cage_tag_clusters" > "$cage_tag_clusters_bed"

# Extract TSS from annotation
awk '$3 == "mRNA"' "$reference_annotation" | \
awk 'BEGIN{OFS="\t"} {
    gsub(/"/, "", $1);
    split($9, attr, ";");
    for (i in attr) if (attr[i] ~ /^ID=/) id=substr(attr[i], 4);
    if ($7 == "+") print $1, $4-1, $4, id, ".", $7;
    else if ($7 == "-") print $1, $5-1, $5, id, ".", $7;
}' > "$reference_annotation_tss_bed"

echo "Sorting BED files..."
sort -k1,1 -k2,2n "$cage_tag_dominant_ctss_bed" > "$sorted_dominant_ctss"
sort -k1,1 -k2,2n "$cage_tag_clusters_bed" > "$sorted_tc_regions"
sort -k1,1 -k2,2n "$reference_annotation_tss_bed" > "$sorted_annotation_tss"

# ±10 bp window
awk 'BEGIN{OFS="\t"} {
    start = $2 - 10; if (start < 0) start = 0;
    end = $3 + 10;
    print $1, start, end, $4, $5, $6
}' "$sorted_annotation_tss" > "$reference_annotation_tss_10_bed"

sort -k1,1 -k2,2n "$reference_annotation_tss_10_bed" > "$sorted_annotation_tss_10"

echo "Running bedtools comparisons..."
bedtools intersect -a "$sorted_dominant_ctss" -b "$sorted_annotation_tss" -s -wa -wb > "${directory}exact_dominant_tss_match.bed"
bedtools intersect -a "$sorted_annotation_tss" -b "$sorted_tc_regions" -s -wa -wb > "${directory}tss_overlaps_tc.bed"
bedtools intersect -a "$sorted_dominant_ctss" -b "$sorted_annotation_tss_10" -s -wa -wb > "${directory}fuzzy_dominant_match_10bp.bed"
bedtools closest -a "$sorted_dominant_ctss" -b "$sorted_annotation_tss" -s -D a > "${directory}ctss_to_annotation_distance.bed"

echo ""
echo "Summary:"
echo "--------"
echo "Exact matches (CTSS == TSS): $(wc -l < ${directory}exact_dominant_tss_match.bed)"
echo "Annotated TSS overlapping any TC: $(cut -f1-6 ${directory}tss_overlaps_tc.bed | sort | uniq | wc -l)"
echo "Fuzzy matches (±10bp): $(wc -l < ${directory}fuzzy_dominant_match_10bp.bed)"
echo "Distances file: ${directory}ctss_to_annotation_distance.bed (last column is distance)"
awk -F',' 'NR>1 {sum += $4} END {print "Average TC width:", sum / (NR - 1)}' "$cage_tag_clusters"

# Create gene-level match table
awk 'BEGIN{OFS=","} {
    chrom=$7; tss=$8 + 1; gene=$9; strand=$12;
    sub(/\..*$/, "", gene);
    print gene, tss, strand, chrom
}' "${directory}exact_dominant_tss_match.bed" | sort -u > "${directory}dominant_ctss_exact_tss_table.csv"

echo "Table of dominant CTSSs exactly matching TSS written to: dominant_ctss_exact_tss_table.csv"

# Step 1: Extract promoter-tagged regions with gene info
awk -F',' 'NR>1 && $12 ~ /promoter/ && $13 != "" {
    gsub(/"/, "", $1); gsub(/"/, "", $2); gsub(/"/, "", $3); gsub(/"/, "", $5); gsub(/"/, "", $13);
    OFS="\t"; print $1, $2, $3, $13, ".", $5
}' "$cage_tag_clusters" > "${directory}promoter_tc_regions_raw.bed"
awk '{
    split($4, genes, ";");
    for (g in genes) {
        gene = genes[g];
        sub(/\..*$/, "", gene);
        print gene;
    }
}' "${directory}promoter_tc_regions_raw.bed" | sort -u > "${directory}all_promoter_genes.txt"
bedtools intersect -a "$reference_annotation_tss_bed" -b "${directory}promoter_tc_regions_raw.bed" -s -wa -wb > "${directory}promoter_tss_inside_tc.bed"
cut -f4 "${directory}promoter_tss_inside_tc.bed" | sed 's/\..*//' | sort -u > "${directory}promoter_genes_with_TSS_overlap.txt"
total=$(wc -l < "${directory}all_promoter_genes.txt")
matched=$(wc -l < "${directory}promoter_genes_with_TSS_overlap.txt")
echo "Total expressed genes (promoter-tagged): $total"
echo "Genes with TSS inside their TC: $matched"
if [ "$total" -gt 0 ]; then
  proportion=$(echo "scale=2; $matched / $total * 100" | bc)
  echo "Proportion of expressed genes with TSS support: $proportion%"
else
  echo "No promoter-tagged genes found."
fi
# comm -23 shows lines unique to file 1
comm -23 "${directory}all_promoter_genes.txt" "${directory}promoter_genes_with_TSS_overlap.txt" > "${directory}promoter_genes_without_TSS_overlap.txt"
echo "Missing gene list written to: promoter_genes_without_TSS_overlap.txt"


# Extract genes with exact CTSS == TSS match (gene-level)
cut -f9 "${directory}exact_dominant_tss_match.bed" | sed 's/\..*//' | sort -u > "${directory}genes_with_exact_ctss_match.txt"
total=$(wc -l < "${directory}all_promoter_genes.txt")
exact_match=$(wc -l < "${directory}genes_with_exact_ctss_match.txt")
echo ""
echo "Exact dominant CTSS == TSS matches:"
echo "Total expressed genes (promoter-tagged): $total"
echo "Genes with exact CTSS–TSS match: $exact_match"

if [ "$total" -gt 0 ]; then
  proportion_exact=$(echo "scale=2; $exact_match / $total * 100" | bc)
  echo "Proportion of expressed genes with exact CTSS matching TSS: $proportion_exact%"
else
  echo "No promoter-tagged genes found."
fi
