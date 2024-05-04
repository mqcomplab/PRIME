from modules.esim import *
import numpy as np
import re
import json
import glob

class FrameSimilarity:
    """A class to calculate the similarity between clusters.
    
    Attributes
    ----------
    c0 : numpy.ndarray
        The data for the top cluster.
    input_files : list
        A list of the cluster files.
    summary_file : str
        The path to the summary file.
    n_clusters : int
        The number of clusters to analyze.
    weighted_by_frames : bool
        Whether to weight similarity values by the number of frames.
    n_ary : str
        The similarity metric to use for comparing clusters.
    weight : str
        The weighting scheme to use for comparing clusters.
    sims : dict
        A dictionary to store the similarity values.
    
    Methods
    -------
    calculate_pairwise()
        Calculates pairwise similarity between each cluster and all other clusters.
    calculate_union()
        Calculates the extended similarity between the union of frame in c0 and cluster k.
    calculate_medoid()
        Calculates the pairwise similarity between every frame in c0 and the medoid of each cluster.
    calculate_outlier()
        Calculates the pairwise similarity between every frame in c0 and the outlier of each cluster.
    """
    
    def __init__(self, cluster_folder=None, summary_file=None, trim_frac=None, n_clusters=None, 
                 weighted_by_frames=True, n_ary='RR', weight='nw'):
        """Initializes instances for the FrameSimilarity class.
        
        Parameters
        ----------
        cluster_folder : str
            The path to the folder containing the normalized cluster files.
        summary_file : str
            The path to the summary file containing the number of frames for each cluster.
        trim_frac : float
            The fraction of outliers to trim from the top cluster.
        n_clusters : int
            The number of clusters to analyze, None for all clusters.
        weighted_by_frames : bool
            Whether to weight similarity values by the number of frames.
        n_ary : str
            The similarity metric to use for comparing clusters.
        weight : str
            The weighting scheme to use for comparing clusters.
        
        Notes
        -----
        - For each cluster file, loads the data and calculates the similarity score 
        with the top (c0) cluster.
        - The esim index used is defined by the `n_ary` parameter.
        """
        self.c0 = np.load(f"{cluster_folder}/normed_clusttraj.c0.npy")
        if trim_frac:
            self.c0 = trim_outliers(self.c0, trim_frac=trim_frac, n_ary=n_ary, weight=weight)
        self.input_files = sorted(glob.glob(f"{cluster_folder}/normed_clusttraj.c*"), 
                                  key=lambda x: int(re.findall("\d+", x)[0]))[1:]
        self.summary_file = summary_file
        self.n_clusters = n_clusters
        self.weighted_by_frames = weighted_by_frames
        self.n_ary = n_ary
        self.weight = weight
        self.sims = {}
    
    def calculate_pairwise(self):
        """The similarity score is calculated as the average of pairwise similarity 
        values between each frame in the cluster and the top c0 cluster.
        
        Returns
        -------
        If `frame_weighted_sim` returns `False`, 
            nw_dict (dict): unweighted average similarity values.
        If `frame_weighted_sim` returns `True`,
            w_dict (dict): calls `weight_dict` function to weight similarity values.
        """
        for each, file in enumerate(self.input_files):
            ck = np.load(file)
            self.sims[each] = {}
            for i, x in enumerate(self.c0):
                total = 0
                for j, y in enumerate(ck):
                    c_total = np.sum(np.array([x, y]), axis=0)
                    pair_sim = SimilarityIndex(c_total, 2, n_ary=self.n_ary, weight=self.weight, 
                                               c_threshold=None, w_factor="fraction")()
                    total += pair_sim
                avg = total / len(ck)
                if f"f{i}" not in self.sims[each]:
                    self.sims[each][f"f{i}"] = []
                self.sims[each][f"f{i}"] = avg

        nw_dict = _format_dict(self.sims)
        if not self.weighted_by_frames:
            return nw_dict
        elif self.weighted_by_frames:
            return weight_dict(file_path=None, summary_file=self.summary_file, dict=nw_dict, 
                               n_clusters=self.n_clusters)

    def calculate_union(self):
        """The similarity score is calculated as the union similarity between 
        all frames in the cluster and the top c0 cluster.
        
        Returns
        -------
        If `frame_weighted_sim` returns `False`, 
            nw_dict (dict): unweighted average similarity values.
        If `frame_weighted_sim` returns `True`,
            w_dict (dict): calls `weight_dict` function to weight similarity values.
        """
        for each, file in enumerate(self.input_files):
            ck = np.load(file)
            self.sims[each] = {}
            for i, x in enumerate(self.c0):
                c_total = np.sum(ck, axis=0) + x
                n_fingerprints = len(ck) + 1
                index = SimilarityIndex(c_total, n_fingerprints, n_ary=self.n_ary, weight=self.weight,
                                        c_threshold=None, w_factor="fraction")()
                if f"f{i}" not in self.sims[each]:
                    self.sims[each][f"f{i}"] = []
                self.sims[each][f"f{i}"] = index
        
        nw_dict = _format_dict(self.sims)
        if not self.weighted_by_frames:
            return nw_dict
        elif self.weighted_by_frames:
            return weight_dict(file_path=None, summary_file=self.summary_file, dict=nw_dict, 
                               n_clusters=self.n_clusters)

    def _perform_calculation(self, index_func):
        """Auxillary function for `calculate_medoid` and `calculate_outlier`.

        Parameters
        ----------
        index_func : function
            The function to calculate the medoid or outlier of each cluster.
        
        Returns
        -------
        nw_dict (dict): unweighted average similarity values.
        """
        for each, file in enumerate(self.input_files):
            ck = np.load(file)
            index = index_func(ck, n_ary=self.n_ary, weight=self.weight)
            medoid = ck[index]
            self.sims[each] = {}
            for i, x in enumerate(self.c0):
                c_total = medoid + x
                if f"f{i}" not in self.sims[each]:
                    self.sims[each][f"f{i}"] = []
                pair_sim = SimilarityIndex(c_total, 2, n_ary=self.n_ary, weight=self.weight,
                                           c_threshold=None, w_factor="fraction")()
                self.sims[each][f"f{i}"] = pair_sim
        nw_dict = _format_dict(self.sims)
        return nw_dict
    
    def calculate_medoid(self):
        """The pairwise similarity value between each frame in c0 and the medoid of each cluster 
        is calculated using similarity indices. Calls the `_perform_calculation` aux function.
        
        Returns
        -------
        If `frame_weighted_sim` returns `False`, 
            nw_dict (dict): unweighted average similarity values.
        If `frame_weighted_sim` returns `True`,
            w_dict (dict): calls `weight_dict` function to weight similarity values.
        """
        nw_dict = self._perform_calculation(calculate_medoid)
        if not self.weighted_by_frames:
            return nw_dict
        elif self.weighted_by_frames:
            return weight_dict(file_path=None, summary_file=self.summary_file, dict=nw_dict, 
                               n_clusters=self.n_clusters)
        
    def calculate_outlier(self):
        """The pairwise similarity value between each frame in c0 and the outlier of each cluster 
        is calculated using similarity indices. Calls the `_perform_calculation` auxillary function.
        
        Returns
        -------
        If `frame_weighted_sim` returns `False`, 
            nw_dict (dict): unweighted average similarity values.
        If `frame_weighted_sim` returns `True`,
            w_dict (dict): calls `weight_dict` function to weight similarity values.
        """
        nw_dict = self._perform_calculation(calculate_outlier)
        if not self.weighted_by_frames:
            return nw_dict
        elif self.weighted_by_frames:
            return weight_dict(file_path=None, summary_file=self.summary_file, dict=nw_dict, 
                               n_clusters=self.n_clusters)

def trim_outliers(total_data, trim_frac=0.1, n_ary='RR', weight='nw', removal='nan'):
    """Trims a desired percentage of outliers (most dissimilar) from the dataset 
    by calculating largest complement similarity.

    Parameters
    ----------
    total_data : numpy.ndarray
        The data to trim outliers from.
    trim_frac : float, optional
        The fraction of outliers to trim. The default is 0.1.
    n_ary : str, optional
        The n-ary method. The default is 'RR'.
    weight : str, optional
        The weight method. The default is 'nw'.
    removal : str, optional
        The method of removal. The default is 'nan'.

    Returns
    -------
    numpy.ndarray
        The trimmed dataset.
    """
    n_fingerprints = len(total_data)
    c_total = np.sum(total_data, axis = 0)
    comp_sims = []
    for i, pixel in enumerate(total_data):
        c_total_i = c_total - total_data[i]
        sim_index = SimilarityIndex(c_total_i, n_fingerprints - 1, n_ary=n_ary, weight=weight,
                                    c_threshold=None, w_factor="fraction")()
        comp_sims.append(sim_index)
    comp_sims = np.array(comp_sims)
    cutoff = int(np.floor(n_fingerprints * float(trim_frac)))
    highest_indices = np.argpartition(-comp_sims, cutoff)[:cutoff]
    if removal == 'nan':
        total_data[highest_indices] = np.nan
    elif removal == 'delete':
        total_data = np.delete(total_data, highest_indices, axis=0)
    return total_data

def weight_dict(file_path=None, summary_file=None, dict=None, n_clusters=None):
    """Calculates frame-weighted similarity values by the number of frames in each cluster.

    Parameters
    ----------
    file_path : str, optional
        The path to the file containing the similarity values. The default is None.
    summary_file : str
        The path to the summary file containing the number of frames for each cluster.
    dict : dict
        A dictionary containing the similarity values. The default is None.
    n_clusters : int
        The number of clusters to analyze. The default is None.

    Returns
    -------
    dict
        A dictionary containing the weighted similarity values.
    """
    if file_path:
        with open(file_path, 'r') as file:
            dict = json.load(file)
    elif dict:
        dict = dict
    for key in dict:
        dict[key].pop()
    
    num = np.loadtxt(summary_file, unpack=True, usecols=(1), skiprows=(1), delimiter=',')
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

def _format_dict(dict):
    """Sorts dict to have frame # as the key and attaches the average value to the end of each key.

    Parameters
    ----------
    dict : dict
        A dictionary containing the similarity values.
    
    Returns
    -------
    dict
        A dictionary with the frame number as the key and the average value attached to the end.
    """
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