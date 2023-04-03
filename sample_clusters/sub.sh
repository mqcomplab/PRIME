#!/bin/bash
#SBATCH --job-name=normalize
#SBATCH --output=py.out
#SBATCH --error=py.err
#SBATCH --mail-type=NONE
#SBATCH --mail-user=some_user@some_domain.com
#SBATCH --time=7-00:00:00 # adjust time
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1    # number of CPU core to use
#SBATCH --mem-per-cpu=100gb
#SBATCH --account=group
#SBATCH --qos=group

ml python

python normalize.py
