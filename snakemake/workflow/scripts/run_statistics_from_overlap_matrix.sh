#!/bin/bash 
#SBATCH -p ei-short
#SBATCH -o statistics.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk

source ~/.bashrc 
mamba activate /hpc-home/mahony/miniforge3
conda activate miniconda_dna 

python statistics_from_overlap_matrix.py  
