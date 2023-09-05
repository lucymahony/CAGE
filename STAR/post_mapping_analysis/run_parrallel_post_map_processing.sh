#!/bin/bash -e
#SBATCH -p ei-medium                                              # partition
#SBATCH -n 24                                               # 24 parallel tasks
#SBATCH -c 16                                                # 1 CPU per task (single threaded process)
#SBATCH --mem 102399
#SBATCH -t 0-00:45                                          # maximum of 2 hours
#SBATCH --mail-type=ALL                                     # email me everything
#SBATCH --mail-user=mahony@nbi.ac.uk                    # the email
#SBATCH -o parallel.out   # show me the output
#SBATCH -e parallel.err   # report errors

./list_parrallel_tasks.sh
echo "tissue.list created"
srun post_map_processing.sh


