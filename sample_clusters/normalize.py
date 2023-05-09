"""This script normalizes the data using the min-max normalization method.

Input:
----------------
Break line (-b) - where the line breaks in the input CRD file.
Norm type (-n) - which normalization method to use.
    Variant 2 (v2) optimized for 1-similarity.
    Variant 3 (v3) is normal min-max normalization without optimization.

How it works:
----------------
1. Combine all the files into one files,
2. Normalize all data using the min-max normalization method.
3. Get the min and max for the whole dataset.
4. Concatenate all the normalized data into one file `normed_data.txt`.
5. Normalize each file using the min and max determined from above `normed_clusttraj.c*`.

Example usage:
----------------
>>> python normalize.py -b 65 -n v2
>>> python normalize.py -b 65 -n v3
"""
import sys
sys.path.insert(0, '../')
import modules as mod
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-b', '--break_line', help='break line', required=True, type=int)
parser.add_argument('-n', '--norm_type', help='type of normalization', required=True, type=str)
args = parser.parse_args()

data = mod.read_cpptraj(args.break_line)
norm = mod.Normalize(data=data)
if args.norm_type == 'v2':
    normed_data = norm.get_v2_norm()
elif args.norm_type == 'v3':
    normed_data = norm.get_v3_norm()
min, max = norm.get_min_max()
np.savetxt("normed_data.txt", normed_data)

data = mod.read_cpptraj(args.break_line, min=min, max=max, normalize=True)