#!/bin/bash 
#SBATCH -p ei-medium
#SBATCH -o post_clustering.out
#SBATCH -c 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk
#SBATCH --mem=16G # 256
#SBATCH --time=0-01:12:00


singularity exec feb_252.img Rscript post_clustering.R

