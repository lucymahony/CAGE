#!/bin/bash 
#SBATCH -p ei-medium
#SBATCH -o R_TC.out
#SBATCH -c 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk
#SBATCH --mem=256G #56G # 256
#SBATCH --time=2-00:21:00

singularity exec feb_252.img Rscript generate_TC_bed.R distclu 20 5 no_merge 16
