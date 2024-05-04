"This script aims to find the representative frame for each method below."""
from modules.esim import calculate_medoid
from modules.sim_calc import trim_outliers
import numpy as np
import json
import re

def calculate_max_key(dict):
    """Calculate the key with the maximum value in a dictionary.

    Parameters
    ----------
    dict : dict
        The dictionary to be searched.
    
    Returns
    -------
    max_key : int
        The key with the maximum value.
    """
    max_val = float('-inf')
    max_key = None
    for k, v in dict.items():
        for num in v:
            if num > max_val:
                max_val = num
                max_key = k
    max_key = int(re.findall(r'\d+', max_key)[0])
    return max_key

def gen_all_methods_max(sim_folder='nw', norm_folder='v3_norm', weighted_by_frames=True, trim_frac=0.1, n_ary='RR', weight='nw', output_name='rep'):
    """Generate the representative frame for each method.

    Parameters
    ----------
    weighted_by_frames : bool, optional
        Similarity is weighted by frames. The default is True.
    trim_frac : float, optional
        The fraction of outliers to be trimmed. The default is 0.1.
    n_ary : str, optional
        The n-ary method. The default is 'RR'.
    weight : str, optional
        The weight method. The default is 'nw'.
    output_name : str, optional
        The output name. The default is 'rep'.
    """
    if weighted_by_frames is True:
        w = "w_"
    elif weighted_by_frames is False:
        w = ""
    if trim_frac:
        t = f"_t{int(float(trim_frac) * 100)}"
    elif not trim_frac:
        t= ""
    with open(f"{sim_folder}/{w}{output_name}_{n_ary}{t}.txt","w") as output:
        output.write("# Frame number with max values by method: medoid_all, medoid_c0, medoid_c0(trimmed), pairwise, union, medoid, outlier\n")
        
        # medoid_all
        c_all = np.load(f"{norm_folder}/normed_data.npy")
        output.write(f"{calculate_medoid(c_all, n_ary=n_ary, weight=weight)}, ")
        
        # medoid_c0 (untrimmed)
        c0 = np.load(f"{norm_folder}/normed_clusttraj.c0.npy")
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
        with open(f"{sim_folder}/{w}pairwise_{n_ary}{t}.txt", "r") as file:
            pairwise = json.load(file)
        output.write(f"{calculate_max_key(pairwise)}, ")

        # union
        with open(f"{sim_folder}/{w}union_{n_ary}{t}.txt", "r") as file:
            union = json.load(file)
        output.write(f"{calculate_max_key(union)}, ")

        # medoid
        with open(f"{sim_folder}/{w}medoid_{n_ary}{t}.txt", "r") as file:
            medoid = json.load(file)
        output.write(f"{calculate_max_key(medoid)}, ")
        
        # outlier
        with open(f"{sim_folder}/{w}outlier_{n_ary}{t}.txt", "r") as file:
            outlier = json.load(file)
        output.write(f"{calculate_max_key(outlier)}")

def gen_one_method_max(method, sim_folder='nw', norm_folder='v3_norm', weighted_by_frames=True, trim_frac=0.1, n_ary='RR', weight='nw', output_name='rep'):
    """Generate the representative frame for each method.

    Parameters
    ----------
    method : str
        The method to find the representative frame.
    weighted_by_frames : bool, optional
        Similarity is weighted by frames. The default is True.
    trim_frac : float, optional
        The fraction of outliers to be trimmed. The default is 0.1.
    n_ary : str, optional
        The n-ary method. The default is 'RR'.
    weight : str, optional
        The weight method. The default is 'nw'.
    output_name : str, optional
        The output name. The default is 'rep'.
    
    Raises
    ------
    ValueError
        Invalid method. Choose from 'medoid_all', 'medoid_c0', 
        'medoid_c0(trimmed)', 'pairwise', 'union', 'medoid', 'outlier'
    """
    if method not in ['medoid_all', 'medoid_c0', 'medoid_c0(trimmed)', 'pairwise', 'union', 'medoid', 'outlier']:
        raise ValueError("Invalid method. Choose from 'medoid_all', 'medoid_c0', 'medoid_c0(trimmed)', 'pairwise', 'union', 'medoid', 'outlier'")
    if weighted_by_frames is True:
        w = "w_"
    elif weighted_by_frames is False:
        w = ""
    if trim_frac:
        t = f"_t{int(float(trim_frac) * 100)}"
    elif not trim_frac:
        t= ""
    
    with open(f"{sim_folder}/{w}{output_name}_{n_ary}{t}_{method}.txt","w") as output:
        output.write(f"# Frame number with max values by method: {method}\n")
        
        if method == 'medoid_all':
            c_all = np.load(f"{norm_folder}/normed_data.npy")
            output.write(f"{calculate_medoid(c_all, n_ary=n_ary, weight=weight)}")
        
        elif method == 'medoid_c0':
            c0 = np.load(f"{norm_folder}/normed_clusttraj.c0")
            output.write(f"{calculate_medoid(c0, n_ary=n_ary, weight=weight)}")
        
        elif method == 'medoid_c0(trimmed)':
            if not trim_frac:
                raise ValueError("No trimmed frac available")
            elif trim_frac:
                trim_c0 = trim_outliers(c0, trim_frac=trim_frac, n_ary=n_ary, weight=weight, removal='delete')
                index = calculate_medoid(trim_c0)
                search = trim_c0[index]
                new_index = np.where((c0 == search).all(axis=1))[0]
                output.write(f"{new_index[0]}")
        
        prime_methods = ['pairwise', 'union', 'medoid', 'outlier']

        if method in prime_methods:
            with open(f"{sim_folder}/{w}{method}_{n_ary}{t}.txt", "r") as file:
                data = json.load(file)
            output.write(f"{calculate_max_key(data)}")