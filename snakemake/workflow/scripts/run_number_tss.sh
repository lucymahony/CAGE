#!/bin/bash 
#SBATCH -p ei-medium
#SBATCH -o number_tss.out
#SBATCH -c 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk
#SBATCH --mem=10G
#SBATCH --time=0-0:12:00

my_annotation_file='/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/tmp/snakemake_intermediate_data/no_merge_reps/finalised_custom_annotations/ALL_TC_merged_output.csv'
output_csv_file='/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/tmp/snakemake_intermediate_data/finalised_custom_annotations_ALL_TCs.csv'
output_directory='/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/tmp/snakemake_intermediate_data/'


source ~/.bashrc 
mamba activate /hpc-home/mahony/miniforge3
conda activate miniconda_dna 
python number_tss_to_expression_level.py ${my_annotation_file} ${output_csv_file} ${output_directory}

