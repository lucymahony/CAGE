#!/bin/bash 
#SBATCH -p ei-medium
#SBATCH -o count_enhancers.out
#SBATCH -c 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk
#SBATCH --mem=10G
#SBATCH --time=0-01:12:00



# Input files
source package b0ed0698-358b-4c9b-9d21-603ea8d6e478  # bedtools 2.31.0
file_directory="/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data"

#Array of BED files
enhancer_files=("shared_L_cadenza_at_least_2_overlaps_max_score_enhancers.bed" "shared_R_cadenza_at_least_2_overlaps_max_score_enhancers.bed" "shared_RO_fielder_at_least_2_overlaps_max_score_enhancers.bed" "shared_LE_fielder_at_least_2_overlaps_max_score_enhancers.bed" "shared_SP_fielder_at_least_2_overlaps_max_score_enhancers.bed" "shared_IS_fielder_at_least_2_overlaps_max_score_enhancers.bed")


# Use bedtools intersect to count how many enhancers overlap in all 6 files, how many dont overlap in any file, and the total number of enhancers in at least one file where overlaps are counted as one enhancer
genes_in_all_six_files=$(bedtools intersect -a "$file_directory/${enhancer_files[0]}" -b "$file_directory/${enhancer_files[1]}" "$file_directory/${enhancer_files[2]}" "$file_directory/${enhancer_files[3]}" "$file_directory/${enhancer_files[4]}" "$file_directory/${enhancer_files[5]}" -wa -u | wc -l)
genes_in_just_one_file=$(bedtools intersect -a "$file_directory/${enhancer_files[0]}" -b "$file_directory/${enhancer_files[1]}" "$file_directory/${enhancer_files[2]}" "$file_directory/${enhancer_files[3]}" "$file_directory/${enhancer_files[4]}" "$file_directory/${enhancer_files[5]}" -wa -v | wc -l)
genesat_least_one_file=$(bedtools intersect -a "$file_directory/${enhancer_files[0]}" -b "$file_directory/${enhancer_files[1]}" "$file_directory/${enhancer_files[2]}" "$file_directory/${enhancer_files[3]}" "$file_directory/${enhancer_files[4]}" "$file_directory/${enhancer_files[5]}" -wa | wc -l)


# Output summary text file
output_summary_text_file="summary_enhancer_analysis.txt"

echo -e "\nGenes present in all 6 files: $genes_in_all_six_files" >> "$output_summary_text_file"
echo "Genes present in only one file: $genes_in_just_one_file" >> "$output_summary_text_file"
echo "Genes present in at least one file / Total number of genes: $genesat_least_one_file" >> "$output_summary_text_file"


echo "Analysis complete"


