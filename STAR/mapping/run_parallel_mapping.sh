#!/bin/bash -e
#SBATCH -p ei-medium                                              # partition
#SBATCH -n 12                                               # 24 parallel tasks
#SBATCH -c 16                                                # 1 CPU per task (single threaded process)
#SBATCH --mem 102399
#SBATCH -t 0-02:00                                          # maximum of 2 hours
#SBATCH --mail-type=ALL                                     # email me everything
#SBATCH --mail-user=mahony@nbi.ac.uk                    # the email
#SBATCH -o run_parallel.out   # show me the output
#SBATCH -e run_parallel.err   # report errors 

srun star_mapping_reads.sh

