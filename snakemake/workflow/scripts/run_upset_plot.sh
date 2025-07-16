#!/bin/bash 
#SBATCH -p ei-medium
#SBATCH -o upset_plot.out
#SBATCH -c 1
#SBATCH --mem=40G # 256
#SBATCH --time=0-06:12:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk
#SBATCH --array=0-1

source ~/.bashrc 
mamba activate /hpc-home/mahony/miniforge3
conda activate miniconda_dna 


base_dir=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/no_merge_reps/

clustering_name_1=paraclu_md100_single_0.3


# All combinations of base_dir and clustering_name
combinations=(
    "$base_dir $clustering_name_1"
)

config="${combinations[$SLURM_ARRAY_TASK_ID]}"
read base_dir clustering_name <<< "$config"

echo "Running task $SLURM_ARRAY_TASK_ID with: $base_dir $clustering_name"

python upset_plot.py $base_dir $clustering_name

