#!/bin/bash 
#SBATCH -p ei-medium
#SBATCH -o snake.out
#SBATCH -c 16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk
#SBATCH --mem=128G

source ~/.bashrc 
mamba activate /hpc-home/mahony/miniforge3
conda run -n snakemake snakemake -j1 --debug-dag --use-conda ../intermediate_data/snakemake_intermediate_data/CS.star.index -s workflow/rules/align.smk --use-conda --conda-frontend conda --cores 16


#conda run -n snakemake snakemake -j1 --debug-dag --use-conda ../intermediate_data/snakemake_intermediate_data/untrimmed_counts_concatenated.txt -s workflow/rules/trim.smk --use-conda --conda-frontend conda --cores 1 --rerun-incomplete

#source activate snake
#snakemake --unlock
#snakemake -j1 -p --use-conda ../intermediate_data/snakemake_intermediate_data/trimmed_SP1_1P.fastq.gz

#snakemake -j1 -p --use-conda ../intermediate_data/snakemake_intermediate_data/trimmed_SP1_1P.count

#snakemake -j16 -p --use-conda --rerun-incomplete ../intermediate_data/snakemake_intermediate_data/trimmed_counts_concatenated.txt

#snakemake -j16 -p --use-conda ../intermediate_data/snakemake_intermediate_data/star_index_201216_Fielder_pseudomolecules_V1+unanchored_contigs/

#snakemake -j16 -p --use-conda ../intermediate_data/snakemake_intermediate_data/bowtie_index_201216_Fielder_pseudomolecules_V1+unanchored_contigs.1.bt2l

#snakemake -j16 -p --use-conda ../intermediate_data/snakemake_intermediate_data/bwa_index_201216_Fielder_pseudomolecules_V1+unanchored_contigs.amb

#snakemake -j16 -p --use-conda ../intermediate_data/snakemake_intermediate_data/hisat_index_201216_Fielder_pseudomolecules_V1+unanchored_contigs

#snakemake -j16 -p --use-conda ../intermediate_data/snakemake_intermediate_data/all_indexes_report_201216_Fielder_pseudomolecules_V1+unanchored_contigs.txt

