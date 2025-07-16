#!/bin/bash 
#SBATCH -p ei-medium
#SBATCH -o custom_annot_2.out
#SBATCH -c 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk
#SBATCH --mem=20G
#SBATCH --time=1-00:12:00

singularity exec feb_252.img Rscript custom_promoter_annotations.R no_merge_reps no_merge_distclu_md20_single_5_myCAGEset.rds SRO1,SRO2,SRO3,SLE1,SLE2,SLE3,SSP1,SSP2,SSP3,SIS1,SIS2,SIS3,CR1,CR2,CR3,CR4,CR5,CL1,CL2,CL3,CL4,CL5
