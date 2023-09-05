#!/bin/bash -e
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH -c 16 #num cores per task
#SBATCH --mem 300G # memory pool for all cores
#SBATCH -o fielder_indexing.out  # STDOUT
#SBATCH -e fielder_indexing.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=mahony@nbi.ac.uk # send-to address
#SBATCH -p ei-medium

source package /tgac/software/testing/bin/STAR-2.3.0 # using version 2.3.0 as later versions seem to have bugs

# Variables
GENOME_DIR="../../intermediate_data/genome_directory" # Directory for storing output genome data
FASTA="../../input_data/fielder_genome_data/201215_Fielder_pseudomolecules_V1+unanchored_contigs.fasta" # Path to reference genome .fasta file
GTF="../../input_data/fielder_genome_data/fielder.release.gtf" # Path to gene annotations .gtf file
THREADS=16 # Number of threads

# Create genome index with STAR
STAR --runThreadN $THREADS --runMode genomeGenerate --genomeDir $GENOME_DIR --genomeFastaFiles $FASTA --sjdbGTFfile $GTF --sjdbOverhang 100 --limitGenomeGenerateRAM 100111047092

echo "Genome indexing has finished"
