#!/bin/bash
#SBATCH --job-name=3N0957s1
#SBATCH --output=py.out
#SBATCH --error=py2.err
#SBATCH --mail-type=NONE
#SBATCH --mail-user=some_user@some_domain.com
#SBATCH --time=7-00:00:00 # adjust time
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1    # number of CPU core to use
#SBATCH --mem-per-cpu=100gb
#SBATCH --account=rmirandaquintana
#SBATCH --qos=rmirandaquintana

ml python
python scripts/rep_frames.py -s nw -n v3_norm
python scripts/rep_frames.py -s nw -n v3_norm -t 0.1
python scripts/rep_frames.py -s nw -n v3_norm -t 0.2
python scripts/rep_frames.py -s nw -n v3_norm -t 0.3
python scripts/rep_frames.py -s nw -n v3_norm -i SM
python scripts/rep_frames.py -s nw -n v3_norm -t 0.1 -i SM
python scripts/rep_frames.py -s nw -n v3_norm -t 0.2 -i SM
python scripts/rep_frames.py -s nw -n v3_norm -t 0.3 -i SM
