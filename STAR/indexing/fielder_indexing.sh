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

source package 266730e5-6b24-4438-aecb-ab95f1940339 #STAR
source package 638df626-d658-40aa-80e5-14a275b7464b #Samtools 1.15.1

# Variables
GENOME_DIR="" # Directory for storing output genome data
FASTA="../../input_data/fielder_genome_data/201216_Fielder_pseudomolecules_V1+unanchored_contigs.fasta" # Path to reference genome .fasta file
GTF="../../input_data/fielder_genome_data/fielder.release.gtf" # Path to gene annotations .gtf file
THREADS=16 # Number of threads

# Create genome index with STAR
STAR --runThreadN $THREADS --runMode genomeGenerate --genomeDir $GENOME_DIR --genomeFastaFiles $FASTA --sjdbGTFfile $GTF --sjdbOverhang 100 --limitGenomeGenerateRAM 100111047092

echo "Genome indexing has finished"
