#!/bin/bash -e
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH -c 2 #num cores per task
#SBATCH --mem 200000 # memory pool for all cores
#SBATCH -o stats.out  # STDOUT
#SBATCH -e stats.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=mahony@nbi.ac.uk # send-to address
#SBATCH -p ei-medium
#SBATCH -t 0-01:00  

source package 638df626-d658-40aa-80e5-14a275b7464b #Samtools 1.15.1

input_sam="/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/mapped_reads/IS2Aligned.out.sam"
awk '($1 !~ /^@/) {print $5}' "$input_sam" | sort | uniq -c
samtools view "$input_sam" | \
awk 'BEGIN{FS="\t"} {for(i=1; i<=NF; i++) if ($i ~ /NH:i:/) print $i}' | \
cut -d':' -f3 | \
sort | \
uniq -c

