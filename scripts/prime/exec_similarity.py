"""PRIME is a command line tool for structure prediction

-h - for help with the argument options.
-m - methods (required)
-n - number of clusters (required)
-i - similarity index (required)
-t - Fraction of outliers to trim in decimals.
-w - Weighing clusters by frames it contains.
-d - directory where the `normed_clusttraj.c*` files are located.
-s - location where summary with population file is located

Potential Error:
 w_dict[key] = [old_list[i] * v for i, v in enumerate(weights)]
IndexError: list index out of range
Check -n argument to ensure it is not greater than the number of clusters in the directory.
"""
import os

os.system('python ../../utils/similarity.py -m union -n 10 -i SM -t 0.1  \
          -d ../normalization -s ../clusters/outputs/summary_20.txt')