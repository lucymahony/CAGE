#!/bin/bash -e
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH -c 16 #num cores per task
#SBATCH --mem 200000 # memory pool for all cores
#SBATCH -o x.out  # STDOUT
#SBATCH -e x.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=mahony@nbi.ac.uk # send-to address
#SBATCH -p ei-medium


tissues=("IS1" "IS2" "IS3" "SP1" "SP2" "SP3" "RO1" "RO2" "RO3" "LE1" "LE2" "LE3")
types=("uniquelymapped" "Aligned.out")

# Loop over each file type
for type in "${types[@]}"; do
    # Loop over each tissue type
    for tissue in "${tissues[@]}"; do
        # Formulate the SAM file name based on the type and tissue
        if [[ "$type" == "uniquelymapped" ]]; then
            SAM_FILE="${type}_${tissue}.sam"
        else
            SAM_FILE="${tissue}${type}.sam"
        fi
        
        # Check if the file exists before proceeding
        if [ ! -f "$SAM_FILE" ]; then
            echo "File $SAM_FILE not found!"
            continue
        fi

        # Extract the XM:i: values, get the number part, and then create a histogram
        grep "XM:i:" $SAM_FILE | awk -F: '{print $NF}' | sort | uniq -c > "${tissue}_${type}_XM_histogram.txt"

        echo "The distribution of the number of mismatches for $SAM_FILE has been saved to ${tissue}_${type}_XM_histogram.txt"
    done
done

echo "Processing completed."

