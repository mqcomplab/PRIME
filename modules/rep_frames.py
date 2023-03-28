"""This script aims to find the representative frame for each method below. """
from modules.sim_modules_vector import calculate_medoid
from modules.sim_calc import trim_outliers
import numpy as np
import json
import re

def calculate_max_key(dict):
    """ Find the key with the max value """
    max_val = float('-inf')
    max_key = None

    for k, v in dict.items():
        for num in v:
            if num > max_val:
                max_val = num
                max_key = k

    max_key = int(re.findall(r'\d+', max_key)[0])
    return max_key

def gen_method_max(weighted=True, trim=0.1, n_ary="RR", weight='nw', output_name="rep"):
    if weighted is True:
        w = "w_"
    elif weighted is False:
        w = ""
    if trim is not None:
        t = f"_t{int(float(trim) * 100)}"
    elif trim is None:
        t= ""
    with open(f"{w}{output_name}_{n_ary}{t}.txt","w") as output:
        output.write("# Frame number with max values by method: medoid_all, medoid_c0, pairwise, union, medoid, outlier\n")
        
        c_all = np.genfromtxt("new_clusters/normed_data.txt")
        output.write(f"{calculate_medoid(c_all, n_ary=n_ary, weight=weight)}, ")

        c0 = np.genfromtxt("new_clusters/normed_clusttraj.c0")
        if trim is None:
            output.write(f"{calculate_medoid(c0, n_ary=n_ary, weight=weight)}, ")
        if trim is not None:
            trim_c0 = trim_outliers(c0, trim_frac=trim, n_ary=n_ary, weight=weight)
            index = calculate_medoid(trim_c0)
            search = trim_c0[index]
            new_index = np.where((c0 == search).all(axis=1))[0]
            output.write(f"{new_index[0]}, ")
        
        with open(f"{w}pairwise_{n_ary}{t}.txt", "r") as file:
            pairwise = json.load(file)
        output.write(f"{calculate_max_key(pairwise)}, ")

        with open(f"{w}union_{n_ary}{t}.txt", "r") as file:
            union = json.load(file)
        output.write(f"{calculate_max_key(union)}, ")

        with open(f"{w}medoid_{n_ary}{t}.txt", "r") as file:
            medoid = json.load(file)
        output.write(f"{calculate_max_key(medoid)}, ")
        
        with open(f"{w}outlier_{n_ary}{t}.txt", "r") as file:
            outlier = json.load(file)
        output.write(f"{calculate_max_key(outlier)}")

