#!/bin/bash -e
#SBATCH -p ei-medium                                              # partition
#SBATCH -n 1                                               # 24 parallel tasks
#SBATCH -c 1                                                # 1 CPU per task (single threaded process)
#SBATCH --mem 102399
#SBATCH -t 0-03:45                                          # maximum of 2 hours
#SBATCH --mail-type=ALL                                     # email me everything
#SBATCH --mail-user=mahony@nbi.ac.uk                    # the email
#SBATCH -o number_reads_mapping.out   # show me the output
#SBATCH -e number_reads_mapping.err   # report errors

# Source Samtools 1.15.1
source package 638df626-d658-40aa-80e5-14a275b7464b

# Path to SAM files
intermediate_directory="/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/mapped_reads"

# Initialize the CSV file with the header
echo "tissue,number_of_hits_in_reference,number_of_reads" > number_reads_mapping.csv

# List of tissue samples
tissues=("IS1" "IS2" "SP1" "SP2" "SP3" "RO1" "RO2" "RO3" "LE1" "LE2" "LE3")

# Loop over each tissue and process the corresponding SAM file
for tissue in "${tissues[@]}"; do
  input_sam="$intermediate_directory/${tissue}Aligned.out.sam"

  # Check if the SAM file exists
  if [[ -f $input_sam ]]; then
    # Use samtools and awk to extract the data
    samtools view "$input_sam" | \
    awk 'BEGIN{FS="\t"} {for(i=1; i<=NF; i++) if ($i ~ /NH:i:/) print $i}' | \
    cut -d':' -f3 | \
    sort | \
    uniq -c | \
    # Parse the output and append it to the CSV file
    awk -v tissue="$tissue" '{print tissue "," $2 "," $1}' >> number_reads_mapping.csv
  else
    echo "SAM file for $tissue not found"
  fi
done

