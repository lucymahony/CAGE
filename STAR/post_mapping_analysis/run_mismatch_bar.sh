#!/bin/bash -e
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH -c 1 #num cores per task
#SBATCH --mem 200000 # memory pool for all cores
#SBATCH -o misbar.out  # STDOUT
#SBATCH -e misbar.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=mahony@nbi.ac.uk # send-to address
#SBATCH -p ei-medium


source ~/.bashrc
python3 -u bar_chart_mismatches.py

