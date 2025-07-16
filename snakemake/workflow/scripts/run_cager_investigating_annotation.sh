#!/bin/bash 
#SBATCH -p ei-medium
#SBATCH -o cager_annot.out
#SBATCH -c 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk
#SBATCH --mem=100G # 256
#SBATCH --time=1-00:12:00

singularity exec nov_24.img Rscript cager_investigating_annotation.R
