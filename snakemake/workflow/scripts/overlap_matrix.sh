#!/bin/bash
#SBATCH -p ei-medium
#SBATCH -o overlap_matrix.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk

source package b0ed0698-358b-4c9b-9d21-603ea8d6e478  # bedtools 2.31.0

# Directories
INPUT_DIR="/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/no_merge_reps/finalised_custom_annotations"
WORK_DIR="/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/tmp_overlap_matrix"
mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

ALL="$INPUT_DIR/ALL_TC_merged_output.csv"
ALL_BED="ALL.bed"
FILTERED_ALL="ALL_filtered.csv"

# Step 1: Filter ALL and create ID column
head -n 1 "$ALL" | sed 's/^/id,/' > "$FILTERED_ALL"
tail -n +2 "$ALL" | awk -F',' 'BEGIN{OFS=","} $4 > 0.5 {print $1":"$2"-"$3","$0}' >> "$FILTERED_ALL"

# Step 2: Make BED for overlap
tail -n +2 "$FILTERED_ALL" | awk -F',' 'BEGIN{OFS="\t"} {print $2, $3, $4, $1, ".", $6}' > "$ALL_BED"

# Step 3: Prepare header for matrix
echo -n "id" > matrix.tsv
for FILE in "$INPUT_DIR"/*_TC_merged_output.csv; do
    [[ "$FILE" == "$ALL" ]] && continue
    TISSUE=$(basename "$FILE" | cut -d'_' -f1)
    echo -ne "\t$TISSUE" >> matrix.tsv
done
echo "" >> matrix.tsv

# Step 4: Create id list for join key
cut -d',' -f1 "$FILTERED_ALL" | tail -n +2 > id_list.tmp
cp id_list.tmp matrix.tmp

# Step 5: Generate binary columns
for FILE in "$INPUT_DIR"/*_TC_merged_output.csv; do
    [[ "$FILE" == "$ALL" ]] && continue
    TISSUE=$(basename "$FILE" | cut -d'_' -f1)
    BED_FILE="${TISSUE}.bed"

    tail -n +2 "$FILE" | awk -F',' 'BEGIN{OFS="\t"} $4 > 0.5 {print $1, $2, $3, ".", ".", $5}' > "$BED_FILE"
    bedtools intersect -a "$ALL_BED" -b "$BED_FILE" -s -c | awk '{print ($NF>0 ? 1 : 0)}' > "${TISSUE}.col"
    paste matrix.tmp "${TISSUE}.col" > tmp && mv tmp matrix.tmp
done

# Step 6: Combine id list and matrix
paste id_list.tmp matrix.tmp | sed '1s/^/id\t/' >> matrix.tsv

# Step 7: Merge matrix with filtered ALL CSV
# Replace commas with tabs for join, then back to CSV after
tail -n +1 "$FILTERED_ALL" | tr ',' '\t' > ALL_filtered.tsv
join -t $'\t' -1 1 -2 1 <(sort -k1,1 ALL_filtered.tsv) <(sort -k1,1 matrix.tsv) > joined.tsv

# Restore CSV formatting
# Extract column names from filtered ALL
head -n 1 "$FILTERED_ALL" > header.tmp

# Extract tissue names from matrix header
head -n 1 matrix.tsv | cut -f2- | tr '\t' ',' > matrix_header.tmp

# Combine headers and data
paste -d',' header.tmp matrix_header.tmp > full_header.csv
cat full_header.csv <(cat joined.tsv | tr '\t' ',' ) > ALL_with_overlap_matrix.csv

# Clean up temp header files
rm header.tmp matrix_header.tmp full_header.csv
