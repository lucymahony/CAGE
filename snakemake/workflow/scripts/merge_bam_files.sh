#!/bin/sh
#SBATCH -p ei-medium
#SBATCH -c 1	
#SBATCH --mem 16G				# memory pool for all cores
#SBATCH --time=0-08:00:00				# time limit
#SBATCH --output output_%x		# STDOUT and STDERR
#SBATCH --mail-type=END,FAIL			# notifications for job done & fail
#SBATCH --mail-user=lucy.mahony@earlham.ac.uk	# send-to address


source package c92263ec-95e5-43eb-a527-8f1496d56f1a  # Samtools - 1.18

file_path=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/

samtools merge ${file_path}SRO_merged.bam ${file_path}SRO1.CS.se.star.bam ${file_path}SRO2.CS.se.star.bam ${file_path}SRO3.CS.se.star.bam 

samtools merge ${file_path}SLE_merged.bam ${file_path}SLE1.CS.se.star.bam ${file_path}SLE2.CS.se.star.bam ${file_path}SLE3.CS.se.star.bam

samtools merge ${file_path}SIS_merged.bam ${file_path}SIS1.CS.se.star.bam ${file_path}SIS2.CS.se.star.bam ${file_path}SIS3.CS.se.star.bam

samtools merge ${file_path}SSP_merged.bam ${file_path}SSP1.CS.se.star.bam ${file_path}SSP2.CS.se.star.bam ${file_path}SSP3.CS.se.star.bam

samtools merge ${file_path}CR_merged.bam ${file_path}CR1.CS.se.star.bam ${file_path}CR2.CS.se.star.bam ${file_path}CR3.CS.se.star.bam ${file_path}CR4.CS.se.star.bam ${file_path}CR5.CS.se.star.bam 

samtools merge ${file_path}CL_merged.bam ${file_path}CL1.CS.se.star.bam ${file_path}CL2.CS.se.star.bam ${file_path}CL3.CS.se.star.bam ${file_path}CL4.CS.se.star.bam ${file_path}CL5.CS.se.star.bam
