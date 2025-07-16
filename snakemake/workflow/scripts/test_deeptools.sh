#!/bin/bash 
#SBATCH -p ei-medium
#SBATCH -o deep.out
#SBATCH -c 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk
#SBATCH --mem=1G
#SBATCH --time=0-01:12:00

source package /tgac/software/testing/bin/glibc-2.14 
source package 6daf0c37-1c5e-4cd6-9884-2d0b4d5f9d8f 

bamCoverage --version

