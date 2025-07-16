#!/bin/bash 
#SBATCH -p ei-medium
#SBATCH -o width.out
#SBATCH -c 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk
#SBATCH --mem=50G # 256
#SBATCH --time=0-01:12:00

source ~/.bashrc 
mamba activate /hpc-home/mahony/miniforge3
python plot_iqr_widths_distribution.py ../../../intermediate_data/snakemake_intermediate_data/ ../../../intermediate_data/snakemake_intermediate_data/shared_fielder_cadenza_TC_widths
