#!/bin/bash 
#SBATCH -p ei-medium
#SBATCH -o test_sandbox.out
#SBATCH -c 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk
#SBATCH --mem=1G # 256
#SBATCH --time=0-00:01:00

#singularity exec workflow/scripts/cagesoftwareimage.img Rscript workflow/scripts/cage_r_analysis.R ../intermediate_data/snakemake_intermediate_data/CS.star.cage_r_analysis_results.txt


function convert_sci_to_int {
    sci_num=$1
    std_num=$(printf "%.0f" "$sci_num" 2>/dev/null)
    # Check if conversion was successful and the result is an integer
    if [[ "$std_num" =~ ^-?[0-9]+$ ]]; then
        echo "$std_num"
    else
        echo "Error: Failed to convert '$sci_num' to a valid integer" >&2
        exit 1
    fi
}

line='Chr1A	250457	-	1.03688e+06'

# Use a while loop to process each part of the line
updated_line=$(echo "$line" | while read -r col1 col2 col3 col4; do
    # Convert the scientific notation if present
    col4_converted=$(convert_sci_to_int "$col4")
    # Output the updated line
    echo -e "$col1\t$col2\t$col3\t$col4_converted"
done)

echo "$updated_line"
