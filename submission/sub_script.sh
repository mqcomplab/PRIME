#!/bin/bash
#SBATCH --job-name=2N0957s1
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

python scripts/similarity.py -m pairwise -n 11 -i RR -s Cpptraj_linkage_sieve_eps_2/summary
python scripts/similarity.py -m pairwise -n 11 -i RR -t 0.1 -s Cpptraj_linkage_sieve_eps_2/summary
python scripts/similarity.py -m pairwise -n 11 -i RR -t 0.2 -s Cpptraj_linkage_sieve_eps_2/summary
python scripts/similarity.py -m pairwise -n 11 -i RR -t 0.3 -s Cpptraj_linkage_sieve_eps_2/summary
python scripts/similarity.py -m pairwise -n 11 -i SM -s Cpptraj_linkage_sieve_eps_2/summary
python scripts/similarity.py -m pairwise -n 11 -i SM -t 0.1 -s Cpptraj_linkage_sieve_eps_2/summary
python scripts/similarity.py -m pairwise -n 11 -i SM -t 0.2 -s Cpptraj_linkage_sieve_eps_2/summary
python scripts/similarity.py -m pairwise -n 11 -i SM -t 0.3 -s Cpptraj_linkage_sieve_eps_2/summary
