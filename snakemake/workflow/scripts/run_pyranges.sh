#!/bin/bash 
#SBATCH -p ei-medium
#SBATCH -o pyranges.out
#SBATCH -c 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk
#SBATCH --mem=10G # 256
#SBATCH --time=0-00:12:00

source ~/.bashrc 
mamba activate /hpc-home/mahony/miniforge3
conda activate miniconda_dna 
python pyranges_test.py

