#!/bin/bash 
#SBATCH -p ei-short
#SBATCH -o copy_no.out
#SBATCH -c 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk

source ~/.bashrc 
mamba activate /hpc-home/mahony/miniforge3
conda activate miniconda_dna 
python copy_number.py

