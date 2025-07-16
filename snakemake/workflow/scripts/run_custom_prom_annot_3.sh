#!/bin/bash 
#SBATCH -p ei-medium
#SBATCH -o custom_annot_3.out
#SBATCH -c 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk
#SBATCH --mem=20G
#SBATCH --time=1-00:12:00

singularity exec feb_252.img Rscript custom_promoter_annotations.R merge_reps paraclu_md100_single_0.3_myCAGEset.rds SRO,CR
