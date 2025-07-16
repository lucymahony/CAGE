#!/bin/bash 
#SBATCH -p ei-short
#SBATCH -o go_term.out
#SBATCH -c 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk

source ~/.bashrc 
singularity exec r_sing.img Rscript go_term.R


#source ~/.bashrc 
#mamba activate /hpc-home/mahony/miniforge3
#conda activate miniconda_dna 
#python go_term_visualisation.py