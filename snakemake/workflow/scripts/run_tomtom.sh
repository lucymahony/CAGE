#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH -c 1
#SBATCH -p ei-medium 
#SBATCH -o tomtom.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk
#SBATCH --mem=10G

source package  b0b2f14a-a119-41ce-9139-d53136e846c7 # meme 5.1.0  


tomtom -oc tomtom_output \
       -evalue \
       -dist pearson \
       -thresh 0.1 \
       streme_out/streme.txt \
      ../../../intermediate_data/JASPAR_plants.meme




