#!/bin/bash 
#SBATCH -p ei-medium
#SBATCH -o snake.out
#SBATCH -c 16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk
#SBATCH --mem=256G # 256
#SBATCH --time=3-06:42:00

source ~/.bashrc 
mamba activate /hpc-home/mahony/miniforge3
conda run -n snakemake snakemake --use-singularity --cores 16 ../intermediate_data/snakemake_intermediate_data/cadenza.star.all_post_align_summary_statistics.txt

#conda run -n snakemake snakemake --use-singularity --cores 16 ../intermediate_data/snakemake_intermediate_data/SIS1.CS.se.star_plus.unique.bw -s workflow/Snakefile 

#conda run -n snakemake snakemake --dry-run --use-singularity --cores 16 ../intermediate_data/snakemake_intermediate_data/copy_of_CS_se_star_sorted_bw/CS.star.se.CAGE_fight_r_analysis.txt -s workflow/Snakefile 

#conda run -n snakemake snakemake --use-singularity --cores 16  ../intermediate_data/snakemake_intermediate_data/CS.star.cage_r_analysis_results.txt -s workflow/Snakefile 

#conda run -n snakemake snakemake --use-singularity --cores 16 ../intermediate_data/snakemake_intermediate_data/CL1.CS.se.star.unique.bam -s workflow/Snakefile

#conda run -n snakemake snakemake --use-singularity --cores 1 -np ../intermediate_data/snakemake_intermediate_data/cadenza_mappedto_CS.se.star.average_density.bed -s workflow/Snakefile

# conda run -n snakemake snakemake --use-singularity --cores 16  ../intermediate_data/snakemake_intermediate_data/ -s workflow/Snakefile 

