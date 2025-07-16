#!/bin/bash 
#SBATCH -p ei-medium
#SBATCH -o cis_NAT.out
#SBATCH -c 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk
#SBATCH --mem=10G
#SBATCH --time=0-01:12:00

# Load Bedtools
source package b0ed0698-358b-4c9b-9d21-603ea8d6e478  # bedtools 2.31.0


# Input files
file_directory="/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data"
# Gff of high confidence genes
gff_file="/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/data/wheat_genome/CS_v2.1/iwgsc_refseqv2.1_gene_annotation_200916/iwgsc_refseqv2.1_annotation_200916_HC.gff3"
# Tag cluster present in two or more replicates per tissue/genome 
# Define the output summary text file
output_summary_text_file="${file_directory}/cis_nat_summary.txt"

# Convert GFF file to BED format using awk for bedtools compatibility
awk 'BEGIN {OFS="\t"} $3 == "gene" {print $1, $4-1, $5, $9, ".", $7}' "$gff_file" > "${file_directory}/iwgsc_refseqv2.1_HC_genes.bed"

# Extend gene regions upstream by 50 bp at the 5' end
awk 'BEGIN {OFS="\t"} {
    if ($6 == "+") $2 = ($2 - 50 > 0 ? $2 - 50 : 0);
    else if ($6 == "-") $3 += 50;
    print;
}' "${file_directory}/iwgsc_refseqv2.1_HC_genes.bed" > "${file_directory}/iwgsc_refseqv2.1_HC_genes_extended_5prime.bed"

# Array of BED files
tag_cluster_files=("shared_L_cadenza_at_least_2_overlaps_max_score.bed" "shared_R_cadenza_at_least_2_overlaps_max_score.bed" "shared_RO_fielder_at_least_2_overlaps_max_score.bed" "shared_LE_fielder_at_least_2_overlaps_max_score.bed" "shared_SP_fielder_at_least_2_overlaps_max_score.bed" "shared_IS_fielder_at_least_2_overlaps_max_score.bed")

# Initialize a file to store combined gene counts
all_genes_combined="${file_directory}/all_genes_combined_cis_nats.txt"
> "$all_genes_combined"

# Initialize an array to store the overlap counts
overlap_counts=()

# Loop through each tag cluster file
for bed_file in "${tag_cluster_files[@]}"; do
    output_file="${file_directory}/${bed_file%_at_least_2_overlaps_max_score.bed}_cis_nats.bed"
    output_genes_list="${file_directory}/$(basename ${bed_file%_at_least_2_overlaps_max_score.bed}_cis_nats.txt)"

    # Find overlaps using bedtools intersect with opposite strand constraint
    bedtools intersect -a "${file_directory}/iwgsc_refseqv2.1_HC_genes_extended_5prime.bed" -b "${file_directory}/$bed_file"  -s -wa > "$output_file"
    
    # Extract gene ids and count overlaps
    if [[ -s "$output_file" ]]; then
        awk '{print $4}' "$output_file" | sort | uniq > "$output_genes_list"
        overlap_count=$(wc -l < "$output_genes_list")
        echo "Overlap count for $bed_file: $overlap_count"
    else
        overlap_count=0
    fi

    # Append the gene counts to the combined file
    awk '{print $4}' "$output_file" | sort | uniq | awk -v file="$bed_file" '{print file, $1}' >> "$all_genes_combined"
    

    # Store the overlap count
    overlap_counts+=("$overlap_count")
done

# Calculate genes present in all 6 files and only in one
genes_in_all_six_files=$(awk '{print $2}' "$all_genes_combined" | sort | uniq -c | awk '$1 == 6 {print $2}' | wc -l)
genes_in_just_one_file=$(awk '{print $2}' "$all_genes_combined" | sort | uniq -c | awk '$1 == 1 {print $2}' | wc -l)
genes_in_at_least_one_file=$(awk '{print $2}' "$all_genes_combined" | sort | uniq | wc -l)


# Write the results to the summary file
echo -e "\nGenes present in all 6 files: $genes_in_all_six_files" >> "$output_summary_text_file"
echo "Genes present in only one file: $genes_in_just_one_file" >> "$output_summary_text_file"
echo "Genes present in at least one file / Total number of genes: $genes_in_at_least_one_file" >> "$output_summary_text_file"


echo "Analysis complete. Results saved in $output_summary_text_file"




