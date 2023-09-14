#!/bin/bash -e
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH -c 16 #num cores per task
#SBATCH --mem 200000 # memory pool for all cores
#SBATCH -o blast_db.out  # STDOUT
#SBATCH -e blast_db.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=mahony@nbi.ac.uk # send-to address
#SBATCH -p ei-medium

source package 37f0ffda-9f66-4391-87e2-38ccd398861d # blast+ 2.10.0+-4


makeblastdb -dbtype prot -in fielder.release.protein.fa -input_type fasta -blastdb_version 5


