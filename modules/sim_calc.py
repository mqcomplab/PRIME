# change it back to n_clusters=None for the published version. 
from modules.sim_modules_vector import *
import numpy as np
import re
import json
import glob

class SimilarityCalculator:
    def __init__(self, cluster_file_pattern='new_clusters/normed_clusttraj.', summary_file='Cpptraj_linkage_sieve_eps_2/summary', trim_frac=0.1, n_clusters=11, weighted=True, n_ary='RR', weight='nw'):
        """
        cluster_file_pattern is the file path and pattern of the normalized clusters from preprocess.py. summary_file can be found within the CPPTRAJ clustering output.
        trim_frac is the fraction of outliers user wish to trim from the dominant c0 cluster. n_clusters is the number of cluster the user wish to analyze, None will be all clusters. 
        weighted is when similarity values are weighted by the number of frames of the cluster and this is stored in the summary file from clustering.
        n_ary and weight define the similar metric the user wishes to use. All options found in sim_modules_vector under gen_sim_dict.
        """
        self.c0 = np.genfromtxt(f"{cluster_file_pattern}c0")
        if trim_frac:
            self.c0 = trim_outliers(self.c0, trim_frac=trim_frac, n_ary=n_ary, weight=weight)
        self.input_files = sorted(glob.glob(f"{cluster_file_pattern}c*"), key=lambda x: int(re.findall("\d+", x)[0]))[1:]
        self.summary_file = summary_file
        self.n_clusters = n_clusters
        self.weighted = weighted
        self.n_ary = n_ary
        self.weight = weight
        self.sims = {}
    
    def calculate_pairwise(self):
        for each, file in enumerate(self.input_files):
            ck = np.genfromtxt(file)
            self.sims[each] = {}
            for i, x in enumerate(self.c0):
                total = 0
                for j, y in enumerate(ck):
                    c_total = np.sum(np.array([x, y]), axis=0)
                    pair_sim = calculate_counters(c_total, 2, c_threshold=None, w_factor="fraction")
                    if self.n_ary == 'RR':    
                        total += pair_sim["a"] / pair_sim["p"]
                    elif self.n_ary == 'SM':
                        total += (pair_sim["a"] + pair_sim["d"]) / pair_sim["p"]
                avg = total / len(ck)
                if f"f{i}" not in self.sims[each]:
                    self.sims[each][f"f{i}"] = []
                self.sims[each][f"f{i}"] = avg

        nw_dict = sort_dict_add_avg(self.sims)
        if self.weighted is False:
            return nw_dict
        if self.weighted is True:
            return weight_dict(file_path=None, summary_file=self.summary_file, dict=nw_dict, n_clusters=self.n_clusters)

    def calculate_union(self):
        for each, file in enumerate(self.input_files):
            ck = np.genfromtxt(file)
            self.sims[each] = {}
            for i, x in enumerate(self.c0):
                c_total = np.sum(ck, axis=0) + x
                n_fingerprints = len(ck) + 1
                index = gen_sim_dict(c_total, n_fingerprints, c_threshold=None, w_factor="fraction")
                if f"f{i}" not in self.sims[each]:
                    self.sims[each][f"f{i}"] = []
                self.sims[each][f"f{i}"] = index[self.weight][self.n_ary]
        
        nw_dict = sort_dict_add_avg(self.sims)
        if self.weighted is False:
            return nw_dict
        if self.weighted is True:
            return weight_dict(file_path=None, summary_file=self.summary_file, dict=nw_dict, n_clusters=self.n_clusters)

    def calculate_sims(self, index_func):
        for each, file in enumerate(self.input_files):
            ck = np.genfromtxt(file)
            index = index_func(ck, n_ary=self.n_ary, weight=self.weight)
            medoid = ck[index]
            self.sims[each] = {}
            for i, x in enumerate(self.c0):
                c_total = medoid + x
                pair_sim = calculate_counters(c_total, 2, c_threshold=None, w_factor="fraction")
                if f"f{i}" not in self.sims[each]:
                    self.sims[each][f"f{i}"] = []
                if self.n_ary == 'RR':    
                    self.sims[each][f"f{i}"] = pair_sim["a"] / pair_sim["p"]
                elif self.n_ary == 'SM':
                    self.sims[each][f"f{i}"] = (pair_sim["a"] + pair_sim["d"]) / pair_sim["p"]
        nw_dict = sort_dict_add_avg(self.sims)
        return nw_dict
    
    def calculate_medoid(self):
        nw_dict = self.calculate_sims(calculate_medoid)
        if self.weighted is False:
            return nw_dict
        if self.weighted is True:
            return weight_dict(file_path=None, summary_file=self.summary_file, dict=nw_dict, n_clusters=self.n_clusters)
        
    def calculate_outlier(self):
        nw_dict = self.calculate_sims(calculate_outlier)
        if self.weighted is False:
            return nw_dict
        if self.weighted is True:
            return weight_dict(file_path=None, summary_file=self.summary_file, dict=nw_dict, n_clusters=self.n_clusters)

def trim_outliers(total_data, trim_frac=0.1, n_ary='RR', weight='nw'):
    """ Function will trim a desired percentage of outliers from the dataset by calculating largest complement similarity. """
    n_fingerprints = len(total_data)
    c_total = np.sum(total_data, axis = 0)
    comp_sims = []
    for i, pixel in enumerate(total_data):
        c_total_i = c_total - total_data[i]
        Indices = gen_sim_dict(c_total_i, n_fingerprints - 1)
        sim_index = Indices[weight][n_ary]
        comp_sims.append(sim_index)
    comp_sims = np.array(comp_sims)
    cutoff = int(np.floor(n_fingerprints * float(trim_frac)))
    highest_indices = np.argpartition(-comp_sims, cutoff)[:cutoff]
    total_data = np.delete(total_data, highest_indices, axis=0)
    # total_data[highest_indices] = np.nan
    return total_data

def weight_dict(file_path=None, summary_file=None, dict=None, n_clusters=None):
    """ Similarity values are weighted by the number of frames of the cluster and this is stored in the summary file from clustering. """
    if file_path is not None:
        with open(file_path, 'r') as file:
            dict = json.load(file)
    elif dict is not None:
        dict = dict
    for key in dict:
        dict[key].pop()
    
    num = np.loadtxt(summary_file, unpack = True, usecols=(1), skiprows=(1))
    if n_clusters:
        num = num[0:n_clusters]
    w_sum = np.sum(num, axis=0)
    weights = num / w_sum
    weights = weights[1:]

    w_dict = {}
    for key in dict:
        old_list = dict[key]
        w_dict[key] = [old_list[i] * v for i, v in enumerate(weights)]
    for k in w_dict:
        average = sum(w_dict[k]) / len(w_dict[k])
        w_dict[k].append(average)

    return w_dict

def sort_dict_add_avg(dict):
    """ The dictionary is organized in order of frames and the average is attached to the end of each key. """
    nw_dict = {}
    for i in sorted(dict):
        for k, v in dict[i].items():
            if k not in nw_dict:
                nw_dict[k] = [None] * len(dict)
            nw_dict[k][i] = v
    for k in nw_dict:
        average = sum(nw_dict[k]) / len(nw_dict[k])
        nw_dict[k].append(average)
    
    return nw_dict

