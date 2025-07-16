#!/bin/bash 
#SBATCH -p ei-medium
#SBATCH -o plot.out
#SBATCH -c 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk
#SBATCH --mem=10G # 256
#SBATCH --time=0-01:12:00

source ~/.bashrc 
mamba activate /hpc-home/mahony/miniforge3
python plot_number_tags_different_locations.py

