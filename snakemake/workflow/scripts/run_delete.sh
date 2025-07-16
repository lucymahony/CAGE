#!/bin/bash 
#SBATCH -p ei-medium
#SBATCH -o delete.out
#SBATCH -c 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk
#SBATCH --mem=256G #56G # 256
#SBATCH --time=2-00:21:00

singularity exec feb_252.img Rscript annotate_bed_tcs.R 
