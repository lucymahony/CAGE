#!/bin/bash 
#SBATCH -p ei-medium
#SBATCH -o classify.out
#SBATCH -c 1
#SBATCH --mem=1G # 256
#SBATCH --time=0-00:12:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk

source ~/.bashrc 
mamba activate /hpc-home/mahony/miniforge3
conda activate miniconda_dna 


genome_annotation_file_path="/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/input_data/chinese_spring_genome_data/iwgsc_refseqv2.1_annotation_200916_HC.gff3"  # Path to the genome annotation file e.g. iwgsc_refseqv2.1_annotation_200916_HC.gff3 
triad_file_path="/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/W10_Orthogroups.full.tsv"
tss_file_path="/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/no_merge_reps/TC_GR_CR1_no_merge_distclu_md50_single_5_myCAGEset_my_annot.csv"  # Path to the TSS file

for threshold in 1 10 30 50 100; 
#for threshold in 10; 
do 
    python classify_tss_alignment_triads.py ${genome_annotation_file_path} ${triad_file_path} ${tss_file_path} ${threshold}
done