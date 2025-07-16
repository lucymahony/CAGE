#!/bin/bash 
#SBATCH -p ei-medium
#SBATCH -o test_r_script.out
#SBATCH -c 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk
#SBATCH --mem=100G # 256
#SBATCH --time=1-00:12:00

#singularity exec nov_24.img Rscript test_r_script.R
singularity exec nov_24.img Rscript cage_r_analysis.R ".CS.se.star.unique.sorted.ctss.n.bed" "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data" "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/snakemake/config/samples.tsv" 
