#!/bin/bash

source package 638df626-d658-40aa-80e5-14a275b7464b # samtools
source package 6394519c-541f-479c-b064-dd0b912eac04 # bedtools

mapfile -t TISSUES < tissue.list
tissue=${TISSUES[$SLURM_PROCID]}

read_1=../../intermediate_data/trimmed_reads/read_1/trimmed_${tissue}_1P.fastq.gz
read_2=../../intermediate_data/trimmed_reads/read_2/trimmed_${tissue}_2P.fastq.gz

THREADS=16
GENOME_DIR=../../intermediate_data/genome_directory
OUT_DIR=../../intermediate_data/mapped_reads


STAR --runThreadN $THREADS --genomeDir $GENOME_DIR --readFilesIn $read_1 $read_2 --readFilesCommand zcat --outFileNamePrefix $OUT_DIR/$tissue --outSAMtype BAM SortedByCoordinate --limitGenomeGenerateRAM 1001110470922 --limitBAMsortRAM 3000000000

echo "Done mapping for $tissue"



