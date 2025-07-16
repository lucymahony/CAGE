#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH -c 1
#SBATCH -p ei-medium 
#SBATCH -o logo.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk
#SBATCH --mem=10G


source ~/.bashrc
mamba activate /hpc-home/mahony/miniforge3 
python streme_to_logo.py streme_out/streme.txt streme_out/high_res_motifs
