#!/bin/bash 
#SBATCH -p ei-medium
#SBATCH -o generate_ALL_TC.out
#SBATCH -c 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk
#SBATCH --mem=256G
#SBATCH --time=1-00:12:00

singularity exec feb_252.img Rscript generate_ALL_TC.R distclu 20 5 no_merge 1

