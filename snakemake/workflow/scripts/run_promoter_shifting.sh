#!/bin/bash 
#SBATCH -p ei-medium
#SBATCH -o shifting.out
#SBATCH -c 16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk

singularity exec nov_24.img Rscript promoter_shifting.R

