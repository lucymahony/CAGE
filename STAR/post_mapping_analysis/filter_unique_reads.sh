#!/bin/bash -e
#SBATCH -p ei-medium                                             
#SBATCH -n 1                                               
#SBATCH -c 1                                                
#SBATCH --mem 102399
#SBATCH -t 0-08:00                                          
#SBATCH --mail-type=ALL                                     
#SBATCH --mail-user=mahony@nbi.ac.uk                    
#SBATCH -o filter_reads.out  
#SBATCH -e filter_reads.err   

# Source Samtools 1.15.1
source package 638df626-d658-40aa-80e5-14a275b7464b

# Path to SAM files
intermediate_directory="/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/mapped_reads"


# List of tissue samples
tissues=("IS1" "IS2" "SP1" "SP2" "SP3" "RO1" "RO2" "RO3" "LE1" "LE2" "LE3")

# Loop over each tissue and process the corresponding SAM file
for tissue in "${tissues[@]}"; do
    # Mapped_read_file_name
    input_sam="$intermediate_directory/${tissue}Aligned.out.sam"

    egrep '\<NH:i:[1]{1}\>' "${input_sam}" | cat > uniquelymapped_${tissue}.sam

done

