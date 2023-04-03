import sys
sys.path.insert(0, '../')
import modules as mod
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-b', '--break_line', help='break line', required=True, type=int)
args = parser.parse_args()

# Combine all the files into one files,
# Normalize all data using the min-max normalization,
# and getting the min and max for the whole dataset.

data = mod.read_cpptraj(args.break_line)
norm = mod.Normalize(data=data)
normed_data = norm.get_normed_data()
min, max = norm.get_min_max()
np.savetxt("normed_data.txt", normed_data)

# Normalize each file using the min and max determined from above. 
data = mod.read_cpptraj(args.break_line, min=min, max=max, normalize=True)
