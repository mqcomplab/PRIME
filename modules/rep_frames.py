"This script aims to find the representative frame for each method below."""
from modules.esim import calculate_medoid
from modules.sim_calc import trim_outliers
import numpy as np
import json
import re

def calculate_max_key(dict):
    """Find the key with the max value"""
    max_val = float('-inf')
    max_key = None

    for k, v in dict.items():
        for num in v:
            if num > max_val:
                max_val = num
                max_key = k

    max_key = int(re.findall(r'\d+', max_key)[0])
    return max_key

def gen_method_max(weighted_by_frames=True, trim_frac=0.1, n_ary="RR", weight='nw', output_name="rep"):
    """Generate the representative frame for each method.

    Args:
        weighted_by_frames (bool, optional): Similarity is weighted by frames. Defaults to True.
        trim_frac (float, optional): . Defaults to 0.1.
        n_ary (str, optional): _description_. Defaults to "RR".
        weight (str, optional): _description_. Defaults to 'nw'.
        output_name (str, optional): _description_. Defaults to "rep".
    """
    if weighted_by_frames is True:
        w = "w_"
    elif weighted_by_frames is False:
        w = ""
    if trim_frac:
        t = f"_t{int(float(trim_frac) * 100)}"
    elif not trim_frac:
        t= ""
    with open(f"w/{w}{output_name}_{n_ary}{t}.txt","w") as output:
        output.write("# Frame number with max values by method: medoid_all, medoid_c0, medoid_c0(trimmed), pairwise, union, medoid, outlier\n")
        
        # medoid_all
        c_all = np.genfromtxt("new_clusters/normed_data.txt")
        output.write(f"{calculate_medoid(c_all, n_ary=n_ary, weight=weight)}, ")
        
        # medoid_c0 (untrimmed)
        c0 = np.genfromtxt("new_clusters/normed_clusttraj.c0")
        output.write(f"{calculate_medoid(c0, n_ary=n_ary, weight=weight)}, ")
        
        # medoid_c0 (trimmed)
        if not trim_frac:
            output.write(f"{calculate_medoid(c0, n_ary=n_ary, weight=weight)}, ")
        elif trim_frac:
            trim_c0 = trim_outliers(c0, trim_frac=trim_frac, n_ary=n_ary, weight=weight, removal='delete')
            index = calculate_medoid(trim_c0)
            search = trim_c0[index]
            new_index = np.where((c0 == search).all(axis=1))[0]
            output.write(f"{new_index[0]}, ")
        
        # pairwise
        with open(f"w/{w}pairwise_{n_ary}{t}.txt", "r") as file:
            pairwise = json.load(file)
        output.write(f"{calculate_max_key(pairwise)}, ")

        # union
        with open(f"w/{w}union_{n_ary}{t}.txt", "r") as file:
            union = json.load(file)
        output.write(f"{calculate_max_key(union)}, ")

        # medoid
        with open(f"w/{w}medoid_{n_ary}{t}.txt", "r") as file:
            medoid = json.load(file)
        output.write(f"{calculate_max_key(medoid)}, ")
        
        # outlier
        with open(f"w/{w}outlier_{n_ary}{t}.txt", "r") as file:
            outlier = json.load(file)
        output.write(f"{calculate_max_key(outlier)}")
