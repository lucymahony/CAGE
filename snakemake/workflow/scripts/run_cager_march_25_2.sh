#!/bin/bash 
#SBATCH -p ei-long
#SBATCH -o R_2.out
#SBATCH -c 16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk
#SBATCH --mem=256G #56G # 256
#SBATCH --time=3-12:01:00

singularity exec feb_252.img Rscript generate_transcriptional_clusters.R paraclu 500 0.5 no_merge 16
