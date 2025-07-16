#!/bin/bash 
#SBATCH -p ei-medium
#SBATCH -o R.out
#SBATCH -c 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk
#SBATCH --mem=5G # 256
#SBATCH --time=0-01:12:00

singularity exec feb_25.img Rscript plot_cs_against_cadenza.R
