from archive.sim_modules_vector import *
import numpy as np
import re
import json
import glob

class SimilarityCalculator:
    """A class to calculate the similarity between clusters.
    
    Attributes:
        c0 (numpy.ndarray): The dominant cluster.
        input_files (list): The list of cluster files.
        summary_file (str): The path to the summary file.
        n_clusters (int): The number of clusters to analyze.
        frame_weighted_sim (bool): Whether to weight the similarity values by the number of frames in the cluster.
        n_ary (str): The n_ary similarity metric to use.
        weight (str): The weight to use for the similarity metric.
    
    Methods:
        calculate_pairwise: Calculates the similarity between the dominant cluster and all other clusters.
        calculate_union: Calculates the similarity between the dominant cluster and the union of all other clusters.
        calculate_sims: Calculates the similarity between the dominant cluster and the cluster with the highest similarity to the dominant cluster.
        calculate_medoid: Calculates the similarity between the dominant cluster and the cluster with the lowest average distance to the dominant cluster.
        calculate_outliers: Calculates the similarity between the dominant cluster and the cluster with the highest average distance to the dominant cluster.
    """
    
    def __init__(self, cluster_folder=None, summary_file=None, trim_frac=None, n_clusters=None, frame_weighted_sim=True, n_ary='RR', weight='nw'):
        """Initializes a new instance of the SimilarityCalculator class.
        
        Args:
            cluster_folder (str): The path to the folder containing the normalized cluster files (`preprocess.py`).
            summary_file (str): The path to the summary file containing the number of frames for each cluster (CPPTRAJ clusteringoutput).
            trim_frac (float): The fraction of outliers to trim from the dominant cluster, c0.
            n_clusters (int): The number of clusters to analyze, None for all clusters.
            frame_weighted_sim (bool): Whether to weight similarity values by the number of frames.
            n_ary (str): The similarity metric to use for comparing clusters. 
            weight (str): The weighting scheme to use for comparing clusters.

        Returns:
            None.
            
        Notes:
            Options for n_ary and weight under `sim_modules_vector.py`.
        """
        
        self.c0 = np.genfromtxt(f"{cluster_folder}/normed_clusttraj.c0")
        if trim_frac:
            self.c0 = trim_outliers(self.c0, trim_frac=trim_frac, n_ary=n_ary, weight=weight)
        self.input_files = sorted(glob.glob(f"{cluster_folder}/normed_clusttraj.c*"), key=lambda x: int(re.findall("\d+", x)[0]))[1:]
        self.summary_file = summary_file
        self.n_clusters = n_clusters
        self.frame_weighted_sim = frame_weighted_sim
        self.n_ary = n_ary
        self.weight = weight
        self.sims = {}
    
    def calculate_pairwise(self):
        """ Calculates pairwise similarity between each cluster and all other clusters.

        Notes:
            For each cluster file, loads the data and calculates the similarity score with the dominant c0 cluster.
            The similarity score is calculated as the average of pairwise similarity values between each frame in the cluster and the dominant c0 cluster.
            The similarity metric used is defined by the n_ary parameter and can be either 'RR' or 'SM'.
        
        Returns:
            If `frame_weighted_sim` returns `False`, 
                nw_dict (dict): unweighted average similarity values.
            If `frame_weighted_sim` returns `True`, 
                w_dict (dict): calls `weight_dict` function to weight similarity values.
        """
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
        if self.frame_weighted_sim is False:
            return nw_dict
        if self.frame_weighted_sim is True:
            return weight_dict(file_path=None, summary_file=self.summary_file, dict=nw_dict, n_clusters=self.n_clusters)

    def calculate_union(self):
        """ Calculates the extended similarity between the union of frame in c0 and cluster k.

        Notes:
            For each cluster file, loads the data and calculates the extended similarity.
            The similarity score is calculated as the union similarity between all frames in the cluster and the dominant c0 cluster.
            The similarity metric used is defined by the n_ary parameter and can be either 'RR' or 'SM'.
        
        Returns:
            If `frame_weighted_sim` returns `False`, 
                nw_dict (dict): unweighted average similarity values.
            If `frame_weighted_sim` returns `True`, 
                w_dict (dict): calls `weight_dict` function to weight similarity values.
        """
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
        if self.frame_weighted_sim is False:
            return nw_dict
        if self.frame_weighted_sim is True:
            return weight_dict(file_path=None, summary_file=self.summary_file, dict=nw_dict, n_clusters=self.n_clusters)

    def calculate_sims(self, index_func):
        """ Calculates similarity scores between each cluster and the dominant c0 cluster.

        Args:
            index_func (func): a function that takes a 2D numpy array of fingerprints and returns the index of the 
                        medoid frame, based on the similarity indices used.

        Returns:
            nw_dict (dict): A dictionary containing the average similarity between each pair of clusters.
        
        Notes:
            This function should NOT be used directly. It is called by the `calculate_pairwise` and `calculate_union` functions.
        """
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
    
    def calculate_medoid_c0(self):
        """ Calculates the pairwise similarity between every frame in c0 and the medoid of each cluster.

        Notes:
            Calculate the medoid of each cluster using the `calculate_medoid` function from `sim_modules_vector`.
            The pairwise similarity value between each frame in c0 and the medoid of each cluster is calculated 
            using similarity indices.
            Calls the `calculate_sims` function to calculate the similarity values.
        
        Returns:
            If `frame_weighted_sim` returns `False`, 
                nw_dict (dict): unweighted average similarity values.
            If `frame_weighted_sim` returns `True`, 
                w_dict (dict): calls `weight_dict` function to weight similarity values.
        """
        nw_dict = self.calculate_sims(calculate_medoid)
        if self.frame_weighted_sim is False:
            return nw_dict
        if self.frame_weighted_sim is True:
            return weight_dict(file_path=None, summary_file=self.summary_file, dict=nw_dict, n_clusters=self.n_clusters)
        
    def calculate_outlier_c0(self):
        """ Calculates the pairwise similarity between every frame in c0 and the outlier of each cluster.

        Notes:
            Calculate the outlier of each cluster using the `calculate_outlier` function from `sim_modules_vector`.
            The pairwise similarity value between each frame in c0 and the outlier of each cluster is calculated 
            using similarity indices.
            Calls the `calculate_sims` function to calculate the similarity values.
        
        Returns:
            If `frame_weighted_sim` returns `False`, 
                nw_dict (dict): unweighted average similarity values.
            If `frame_weighted_sim` returns `True`, 
                w_dict (dict): calls `weight_dict` function to weight similarity values.
        """
        nw_dict = self.calculate_sims(calculate_outlier)
        if self.frame_weighted_sim is False:
            return nw_dict
        if self.frame_weighted_sim is True:
            return weight_dict(file_path=None, summary_file=self.summary_file, dict=nw_dict, n_clusters=self.n_clusters)

def trim_outliers(total_data, trim_frac=0.1, n_ary='RR', weight='nw'):
    """ Trims a desired percentage of outliers (most dissimilar) from the dataset by calculating largest complement similarity.

    Args:
        total_data (numpy.ndarray): A 2D array, containing the data to be trimmed.
        trim_frac (float): The fraction of outliers to be removed. Must be between 0 and 1. Defaults to 0.1.
        n_ary (str): The similarity metric to be used. Must be either 'RR' or 'SM'. Defaults to 'RR'.
        weight (str): The weight function to be used. Must be either 'nw' or 'fraction'. Defaults to 'nw'.

    Returns:
        total_data (numpy.ndarray): A 2D array, with the same values as `total_data`, except that a fraction
        `trim_frac` of the rows have been replaced with NaNs, corresponding to the rows with the highest complement
        similarity scores.
    """
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
    # total_data = np.delete(total_data, highest_indices, axis=0)
    total_data[highest_indices] = np.nan
    return total_data

def weight_dict(file_path=None, summary_file=None, dict=None, n_clusters=None):
    """ Calculates frame-weighted similarity values by the number of frames in each cluster.

    Args:
        file_path (str): Path to the json file containing the unweighted similarity values between each pair of clusters.
        summary_file (str): Path to the summary file containing the number of frames in each cluster (CPPTRAJ output).
        dict (dict): A dictionary containing the unweighted similarity values between each pair of clusters.
        n_clusters (int): The number of clusters to analyze. Default is `None`, analyze all clusters from summary file.

    Returns:
        w_dict (dict): frame-weighted similarity values between each pair of clusters.
    """
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
    """ Sorts the dictionary by the key and attaches the average value to the end of each key.

    Args:
        dict (dict): A dictionary containing the similarity values.

    Returns:
        sorted_dict (dict): Sorted by the keys with the average value attached to the end of each key.
    """
    sorted_dict = {}
    for i in sorted(dict):
        for k, v in dict[i].items():
            if k not in sorted_dict:
                sorted_dict[k] = [None] * len(dict)
            sorted_dict[k][i] = v
    for k in sorted_dict:
        average = sum(sorted_dict[k]) / len(sorted_dict[k])
        sorted_dict[k].append(average)   
    return sorted_dict