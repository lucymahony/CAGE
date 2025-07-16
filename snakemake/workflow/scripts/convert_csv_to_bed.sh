#!/bin/sh
#SBATCH -p ei-medium
#SBATCH -c 1	
#SBATCH --mem 10G				# memory pool for all cores
#SBATCH --time=0-01:00:00				# time limit
#SBATCH --output %x.out		# STDOUT and STDERR
#SBATCH --mail-type=END,FAIL			# notifications for job done & fail
#SBATCH --mail-user=lucy.mahony@earlham.ac.uk	# send-to address


INPUT_CSV=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/no_merge_reps/TC_GR_CL1_no_merge_distclu_md20_single_5_myCAGEset_my_annot.csv
OUTPUT_BED=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/no_merge_reps/TC_GR_CL1_no_merge_distclu_md20_single_5_myCAGEset_my_annot.bed


# Extract max score (skip header)
MAX_SCORE=$(awk -F',' 'NR>1 {if ($6 > max) max = $6} END {print max}' "$INPUT_CSV")

# Convert to BED
awk -F',' -v OFS='\t' -v max="$MAX_SCORE" '
BEGIN { getline }  # skip header
NR>1 {
  score = ($6 == "" || max == 0) ? 0 : int(1000 * $6 / max)
  name = "TC_" NR
  print $1, $2, $3, name, score, $5
}' "$INPUT_CSV" > "$OUTPUT_BED"

echo "BED file created: $OUTPUT_BED"
