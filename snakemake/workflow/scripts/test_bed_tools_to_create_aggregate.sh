#!/bin/bash 
#SBATCH -p ei-medium
#SBATCH -o bed.out
#SBATCH -c 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk
#SBATCH --mem=10G
#SBATCH --time=0-01:12:00

# Load Bedtools
source package b0ed0698-358b-4c9b-9d21-603ea8d6e478  # bedtools 2.31.0

# Define input directory and file sets with output categories
input_dir="/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data"

declare -A file_sets
file_sets["SP_fielder"]="SSP1_TC.bed SSP2_TC.bed SSP3_TC.bed"
file_sets["RO_fielder"]="SRO1_TC.bed SRO2_TC.bed SRO3_TC.bed"
file_sets["LE_fielder"]="SLE1_TC.bed SLE2_TC.bed SLE3_TC.bed"
file_sets["IS_fielder"]="SIS1_TC.bed SIS2_TC.bed SIS3_TC.bed"
file_sets["R_cadenza"]="CR1_TC.bed CR2_TC.bed CR3_TC.bed CR4_TC.bed CR5_TC.bed"
file_sets["L_cadenza"]="CL1_TC.bed CL2_TC.bed CL3_TC.bed CL4_TC.bed CL5_TC.bed"


file_sets['fielder_cadenza_all']="shared_RO_fielder_at_least_2_overlaps_max_score.bed shared_LE_fielder_at_least_2_overlaps_max_score.bed shared_IS_fielder_at_least_2_overlaps_max_score.bed shared_SP_fielder_at_least_2_overlaps_max_score.bed shared_R_cadenza_at_least_2_overlaps_max_score.bed shared_L_cadenza_at_least_2_overlaps_max_score.bed"


number_shared_files=6  # Minimum number of files that a region must overlap to be kept

# Loop through each category and process its file set
for category in "${!file_sets[@]}"; do
    files=(${file_sets[$category]})
    output_file="${input_dir}/shared_${category}_at_least_${number_shared_files}_overlaps.bed"
    max_score_file="${input_dir}/shared_${category}_at_least_${number_shared_files}_overlaps_max_score.bed"

    # Initialize an empty array to store final files
    final_files=()

    # Separate the strands (+ and -) and process each separately
    for strand in "+" "-"; do
        # Create temporary files for each strand
        temp_files=()
        for file in "${files[@]}"; do
            input_file="${input_dir}/${file}"
            temp_file="${input_file%.bed}_strand_${strand}.bed"
            awk -v strand="$strand" '$6 == strand' "$input_file" > "$temp_file"
            temp_files+=("$temp_file")
        done

        # Concatenate and sort the files for the current strand
        combined_file="${input_dir}/combined_${strand}_${category}.bed"
        cat "${temp_files[@]}" | sort -k1,1 -k2,2n > "$combined_file"

        # Check if combined file is empty
        if [ ! -s "$combined_file" ]; then
            echo "No regions found on strand $strand for category $category. Skipping..."
            rm "${temp_files[@]}" "$combined_file"
            continue
        fi

        # Find shared regions for the current strand
        shared_file="${input_dir}/shared_${category}_strand_${strand}.bed"
        bedtools intersect -a "$combined_file" -b "${temp_files[@]}" -c | \
            awk -v threshold="$number_shared_files" '$NF >= threshold' > "$shared_file"

        # Check if shared file is empty
        if [ ! -s "$shared_file" ]; then
            echo "No shared regions found on strand $strand for category $category. Skipping..."
            rm "${temp_files[@]}" "$combined_file" "$shared_file"
            continue
        fi

        # Merge overlapping regions with a distance of 10bp, keeping max score and metadata
        final_file="${input_dir}/final_${category}_strand_${strand}.bed"
        bedtools merge -i "$shared_file" -d 10 -c 5,6,7,8 -o collapse,distinct,collapse,collapse | \
            awk 'BEGIN{OFS="\t"} {print $1, $2, $3, ".", $4, $5, $6, $7, $8}' > "$final_file"
        final_files+=("$final_file")

        # Clean up temporary files
        rm "${temp_files[@]}" "$combined_file" "$shared_file"
    done

    # Combine the results from both strands if there are any final files
    if [ ${#final_files[@]} -gt 0 ]; then
        cat "${final_files[@]}" | sort -k1,1 -k2,2n > "$output_file"
        echo "Final shared regions for $category saved to $output_file"
    else
        echo "No shared regions found across all files for category $category."
        continue
    fi

    # Clean up final temporary strand files
    rm -f "${final_files[@]}"

    # Further process the output file to create the max score file
    awk 'BEGIN{OFS="\t"} {
        # Split columns 5, 7, and 8 into arrays
        split($5, scores, ",");
        split($7, starts, ",");
        split($8, ends, ",");

        # Find the maximum score and its index
        max_index = 1;
        max_score = scores[1];
        for (i = 2; i <= length(scores); i++) {
            if (scores[i] > max_score) {
                max_score = scores[i];
                max_index = i;
            }
        }

        # Output the row with the selected values
        print $1, $2, $3, $4, max_score, $6, starts[max_index], ends[max_index];
    }' "$output_file" > "$max_score_file"

    echo "Final shared regions with max score for $category saved to $max_score_file"
done
