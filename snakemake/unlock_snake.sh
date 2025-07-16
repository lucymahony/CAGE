#!/bin/bash 
#SBATCH -p ei-medium
#SBATCH --mem  2G
#SBATCH -o unlock.out
#SBATCH -c 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk


source ~/.bashrc 
mamba activate /hpc-home/mahony/miniforge3
source activate snake
snakemake --unlock
#snakemake --cleanup-metadata ../intermediate_data/snakemake_intermediate_data/*.n.bed

