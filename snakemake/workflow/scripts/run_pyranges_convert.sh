#!/bin/bash 
#SBATCH -p ei-medium
#SBATCH -o pyranges_conv.out
#SBATCH -c 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk
#SBATCH --mem=10G # 256
#SBATCH --time=0-01:12:00

source ~/.bashrc 
mamba activate /hpc-home/mahony/miniforge3
conda activate miniconda_dna 
python turn_bam_to_bigwig_pyranges.py

