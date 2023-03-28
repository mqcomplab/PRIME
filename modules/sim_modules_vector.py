import numpy as np
import random
import glob
import pickle
from math import log
from scipy.stats import rankdata
import time



def calculate_counters(c_total, n_fingerprints, c_threshold=None, w_factor="fraction"):
    """Calculate 1-similarity, 0-similarity, and dissimilarity counters

    Arguments
    ---------
    data_sets : np.ndarray
        Array of arrays. Each sub-array contains m + 1 elements,
        with m being the length of the fingerprints. The first
        m elements are the column sums of the matrix of fingerprints.
        The last element is the number of fingerprints.

    c_threshold : {None, 'dissimilar', int}
        Coincidence threshold.
        None : Default, c_threshold = n_fingerprints % 2
        'dissimilar' : c_threshold = ceil(n_fingerprints / 2)
        int : Integer number < n_fingerprints

    w_factor : {"fraction", "power_n"}
        Type of weight function that will be used.
        'fraction' : similarity = d[k]/n
                     dissimilarity = 1 - (d[k] - n_fingerprints % 2)/n_fingerprints
        'power_n' : similarity = n**-(n_fingerprints - d[k])
                    dissimilarity = n**-(d[k] - n_fingerprints % 2)
        other values : similarity = dissimilarity = 1

    Returns
    -------
    counters : dict
        Dictionary with the weighted and non-weighted counters.

    Notes
    -----
    Please, cite the original papers on the n-ary indices:
    https://jcheminf.biomedcentral.com/articles/10.1186/s13321-021-00505-3
    https://jcheminf.biomedcentral.com/articles/10.1186/s13321-021-00504-4
    """
    # Setting matches
    #total_data = np.sum(data_sets, axis=0)
    #n_fingerprints = int(total_data[-1])
    #c_total = total_data[:-1]
    
    #m = len(c_total)
    
    # Assign c_threshold
    if not c_threshold:
        c_threshold = n_fingerprints % 2
    if c_threshold == 'dissimilar':
        c_threshold = ceil(n_fingerprints / 2)
    if c_threshold == 'min':
        c_threshold = n_fingerprints % 2
    if isinstance(c_threshold, int):
        if c_threshold >= n_fingerprints:
            raise ValueError("c_threshold cannot be equal or greater than n_fingerprints.")
        c_threshold = c_threshold
    if 0 < c_threshold < 1:
        c_threshold *= n_fingerprints
    
    # Set w_factor
    if w_factor:
        if "power" in w_factor:
            power = int(w_factor.split("_")[-1])
            def f_s(d):
                return power**-float(n_fingerprints - d)
    
            def f_d(d):
                return power**-float(d - n_fingerprints % 2)
        elif w_factor == "fraction":
            def f_s(d):
                return d/n_fingerprints
    
            def f_d(d):
                return 1 - (d - n_fingerprints % 2)/n_fingerprints
        else:
            def f_s(d):
                return 1
    
            def f_d(d):
                return 1
    else:
        def f_s(d):
            return 1
    
        def f_d(d):
            return 1
    
    # Calculate a, d, b + c
    new = True
    if new:
        #m = len(c_total)
        a_indices = 2 * c_total - n_fingerprints > c_threshold
        d_indices = n_fingerprints - 2 * c_total > c_threshold
        dis_indices = np.abs(2 * c_total - n_fingerprints) <= c_threshold
        
        a = np.sum(a_indices)
        d = np.sum(d_indices)
        total_dis = np.sum(dis_indices)
        
        a_w_array = f_s(2 * c_total[a_indices] - n_fingerprints)
        d_w_array = f_s(abs(2 * c_total[d_indices] - n_fingerprints))
        total_w_dis_array = f_d(abs(2 * c_total[dis_indices] - n_fingerprints))
        
        w_a = np.sum(a_w_array)
        w_d = np.sum(d_w_array)
        total_w_dis = np.sum(total_w_dis_array)
        
    else:
        a = 0
        w_a = 0
        d = 0
        w_d = 0
        total_dis = 0
        total_w_dis = 0
           
        for s in c_total:
            #print(s, 2*s-n_fingerprints)
            if 2 * s - n_fingerprints > c_threshold:
                a += 1
                w_a += f_s(2 * s - n_fingerprints)
                #print('a')
            elif n_fingerprints - 2 * s > c_threshold:
                d += 1
                w_d += f_s(abs(2 * s - n_fingerprints))
                #print('d')
            else:
                total_dis += 1
                total_w_dis += f_d(abs(2 * s - n_fingerprints))
                #print('dis')
    total_sim = a + d
    total_w_sim = w_a + w_d
    p = total_sim + total_dis
    w_p = total_w_sim + total_w_dis
    
    counters = {"a": a, "w_a": w_a, "d": d, "w_d": w_d,
                "total_sim": total_sim, "total_w_sim": total_w_sim,
                "total_dis": total_dis, "total_w_dis": total_w_dis,
                "p": p, "w_p": w_p}
    return counters
    
def gen_sim_dict(c_total, n_fingerprints, c_threshold=None, w_factor="fraction"):
    counters = calculate_counters(c_total, n_fingerprints, c_threshold=c_threshold, w_factor="fraction")
    # Indices
    # AC: Austin-Colwell, BUB: Baroni-Urbani-Buser, CTn: Consoni-Todschini n
    # Fai: Faith, Gle: Gleason, Ja: Jaccard, Ja0: Jaccard 0-variant
    # JT: Jaccard-Tanimoto, RT: Rogers-Tanimoto, RR: Russel-Rao
    # SM: Sokal-Michener, SSn: Sokal-Sneath n

    # Weighted Indices
    ac_w = (2/np.pi) * np.arcsin(np.sqrt(counters['total_w_sim']/
                                         counters['w_p']))
    bub_w = ((counters['w_a'] * counters['w_d'])**0.5 + counters['w_a'])/\
            ((counters['w_a'] * counters['w_d'])**0.5 + counters['w_a'] + counters['total_w_dis'])
    ct1_w = (log(1 + counters['w_a'] + counters['w_d']))/\
            (log(1 + counters['w_p']))
    ct2_w = (log(1 + counters['w_p']) - log(1 + counters['total_w_dis']))/\
            (log(1 + counters['w_p']))
    ct3_w = (log(1 + counters['w_a']))/\
            (log(1 + counters['w_p']))
    ct4_w = (log(1 + counters['w_a']))/\
            (log(1 + counters['w_a'] + counters['total_w_dis']))
    fai_w = (counters['w_a'] + 0.5 * counters['w_d'])/\
            (counters['w_p'])
    gle_w = (2 * counters['w_a'])/\
            (2 * counters['w_a'] + counters['total_w_dis'])
    ja_w = (3 * counters['w_a'])/\
           (3 * counters['w_a'] + counters['total_w_dis'])
    ja0_w = (3 * counters['total_w_sim'])/\
            (3 * counters['total_w_sim'] + counters['total_w_dis'])
    jt_w = (counters['w_a'])/\
           (counters['w_a'] + counters['total_w_dis'])
    rt_w = (counters['total_w_sim'])/\
           (counters['w_p'] + counters['total_w_dis'])
    rr_w = (counters['w_a'])/\
           (counters['w_p'])
    sm_w =(counters['total_w_sim'])/\
          (counters['w_p'])
    ss1_w = (counters['w_a'])/\
            (counters['w_a'] + 2 * counters['total_w_dis'])
    ss2_w = (2 * counters['total_w_sim'])/\
            (counters['w_p'] + counters['total_w_sim'])


    # Non-Weighted Indices
    ac_nw = (2/np.pi) * np.arcsin(np.sqrt(counters['total_w_sim']/
                                          counters['p']))
    bub_nw = ((counters['w_a'] * counters['w_d'])**0.5 + counters['w_a'])/\
             ((counters['a'] * counters['d'])**0.5 + counters['a'] + counters['total_dis'])
    ct1_nw = (log(1 + counters['w_a'] + counters['w_d']))/\
             (log(1 + counters['p']))
    ct2_nw = (log(1 + counters['w_p']) - log(1 + counters['total_w_dis']))/\
             (log(1 + counters['p']))
    ct3_nw = (log(1 + counters['w_a']))/\
             (log(1 + counters['p']))
    ct4_nw = (log(1 + counters['w_a']))/\
             (log(1 + counters['a'] + counters['total_dis']))
    fai_nw = (counters['w_a'] + 0.5 * counters['w_d'])/\
             (counters['p'])
    gle_nw = (2 * counters['w_a'])/\
             (2 * counters['a'] + counters['total_dis'])
    ja_nw = (3 * counters['w_a'])/\
            (3 * counters['a'] + counters['total_dis'])
    ja0_nw = (3 * counters['total_w_sim'])/\
             (3 * counters['total_sim'] + counters['total_dis'])
    jt_nw = (counters['w_a'])/\
            (counters['a'] + counters['total_dis'])
    rt_nw = (counters['total_w_sim'])/\
            (counters['p'] + counters['total_dis'])
    rr_nw = (counters['w_a'])/\
            (counters['p'])
    sm_nw =(counters['total_w_sim'])/\
           (counters['p'])
    ss1_nw = (counters['w_a'])/\
             (counters['a'] + 2 * counters['total_dis'])
    ss2_nw = (2 * counters['total_w_sim'])/\
             (counters['p'] + counters['total_sim'])

    # Dictionary with all the results
    Indices = {'nw': {'AC': ac_nw, 'BUB':bub_nw, 'CT1':ct1_nw, 'CT2':ct2_nw, 'CT3':ct3_nw,
                      'CT4':ct4_nw, 'Fai':fai_nw, 'Gle':gle_nw, 'Ja':ja_nw,
                      'Ja0':ja0_nw, 'JT':jt_nw, 'RT':rt_nw, 'RR':rr_nw,
                      'SM':sm_nw, 'SS1':ss1_nw, 'SS2':ss2_nw},
                'w': {'AC': ac_w, 'BUB':bub_w, 'CT1':ct1_w, 'CT2':ct2_w, 'CT3':ct3_w,
                      'CT4':ct4_w, 'Fai':fai_w, 'Gle':gle_w, 'Ja':ja_w,
                      'Ja0':ja0_w, 'JT':jt_w, 'RT':rt_w, 'RR':rr_w,
                      'SM':sm_w, 'SS1':ss1_w, 'SS2':ss2_w}}
    return Indices

def calculate_medoid(total_data, n_ary = 'RR', weight = 'nw'):
    """Calculate the medoid of a set"""
    index = len(total_data[0]) + 1
    n_fingerprints = len(total_data)
    min_sim = n_fingerprints * (n_fingerprints + 1) * n_fingerprints
    c_total = np.sum(total_data, axis = 0)
    for i, pixel in enumerate(total_data):
        c_total_i = c_total - total_data[i]
        #data_sets = [np.append(i_sum, len(total_data) - 1)]
        Indices = gen_sim_dict(c_total_i, n_fingerprints - 1)
        sim_index = Indices[weight][n_ary]
        if sim_index < min_sim:
            min_sim = sim_index
            index = i
        else:
            pass
    return index

def calculate_comp_sim_old(total_data, total_sum, c_threshold=None, n_ary = 'RR', weight = 'nw'):
    """Calculate the complementary similarity for all elements"""
    n_fingerprints = len(total_data)
    total = {}
    for i, pixel in enumerate(total_data):
        i_sum = total_sum - total_data[i]
        print(i_sum)
        #data_sets = [np.append(i_sum, len(total_data) - 1)]
        Indices = gen_sim_dict(i_sum, n_fingerprints - 1, c_threshold)
        total[i] = Indices[weight]
    return total

def calculate_comp_sim(total_data, total_sum, c_threshold=None, n_ary = 'RR', weight = 'nw'):
    """Calculate the complementary similarity for all elements"""
    n_fingerprints = len(total_data)
    comp_sums = total_sum - total_data
    total = {}
    for i, c_total in enumerate(comp_sums):
        Indices = gen_sim_dict(c_total, n_fingerprints - 1, c_threshold)
        total[i] = Indices[weight]
    #total = {}
    #for i, pixel in enumerate(total_data):
    #    i_sum = total_sum - total_data[i]
    #    #data_sets = [np.append(i_sum, len(total_data) - 1)]
    #    Indices = gen_sim_dict(i_sum, n_fingerprints - 1, c_threshold)
    #    total[i] = Indices[weight]
    return total

def calculate_comp_sim_rr(total_data, total_sum):
    """Calculate the complementary similarity for all elements"""
    n_fingerprints = len(total_data)
    m = len(total_sum)
    comp_sums = total_sum - total_data
    
    c_threshold = n_fingerprints % 2
    
    #def quick_rr(c_total, n_fingerprints, c_threshold):
    #    return np.sum(((1 + np.sign(2*c_total-n_fingerprints - c_threshold))/2) * (2*c_total-n_fingerprints)/n_fingerprints)
        
    def quick_rr (c_total, n_fingerprints, c_threshold):
        a_indices = 2 * c_total - n_fingerprints > c_threshold
        #print(a_indices)
        a_w_array = (2 * c_total[a_indices] - n_fingerprints)/n_fingerprints
        w_a = np.sum(a_w_array)
        return w_a
        
    total = []
    for i, c_total in enumerate(comp_sums):
        w_a = quick_rr(c_total, n_fingerprints - 1, c_threshold)/m
        total.append((i, w_a))
    total = np.array(total, dtype='float32')
    return total
    
def get_new_index_n(total_data, selected_condensed, n, select_from_n, c_threshold=None, n_ary = 'RR', weight = 'nw'):
    """Select a new diverse molecule"""
    n_total = n + 1
    # min value that is guaranteed to be higher than all the comparisons
    min_value = n * (n + 1) * len(total_data[0])
    
    # placeholder index
    index = len(total_data[0]) + 1
    
    # for all indices that have not been selected
    for i in select_from_n:
        # column sum
        c_total = selected_condensed + total_data[i]
        # calculating similarity
        data_sets = [np.append(c_total, n_total)]
        Indices = gen_sim_dict(data_sets, c_threshold=c_threshold)
        sim_index = Indices[weight][n_ary]
        # if the sim of the set is less than the similarity of the previous diverse set, update min_value and index
        if sim_index < min_value:
            index = i
            min_value = sim_index
    return index

def rank_array(array):
    """Rank an array following the SRD convention.

    Arguments
    ---------
    array : {list, np.array}
        Array whose elements will be ranked.

    Returns
    -------
    Ranked array.

    Notes
    -----
    The SRD convention for ranking is as follows:
    1- Elements are ranked in increasing order (from smallest to biggest).
    2- The rank of elements with the same value is the average of the ranks.
    Example:
    _rank_array([0.1, 0.4, 0.2]) -> np.array([2, 3, 1])
    _rank_array([0.6, 0.2, 0.2, 0.1]) -> np.array([4, 2.5, 2.5, 1])
    """
    order = array.argsort()
    non_repeated_ranks = list(order.argsort() + 1)
    repeated_ranks = list(rankdata(array, method="min"))
    ranks = [0]*len(repeated_ranks)
    for i in range(len(non_repeated_ranks)):
        if repeated_ranks.count(repeated_ranks[i]) == 1:
            ranks[i] = non_repeated_ranks[i]
        else:
            indices = [j for j, x in enumerate(repeated_ranks) if x == repeated_ranks[i]]
            value = 0
            for index in indices:
                value += non_repeated_ranks[index]
            value = value/len(indices)
            for index in indices:
                ranks[index] = value
    return np.array(ranks) - 1

def real_var2_pre(total_data):
    """Pre-processing of normalized continuous data"""
    c_total = np.sum(1 - np.abs(total_data - np.mean(total_data, axis=0)),axis=0)
    return c_total

def calculate_outlier(total_data, n_ary = 'RR', weight = 'nw'):
    """ Calculate the outlier of a set """
    index = len(total_data[0]) + 1
    n_fingerprints = len(total_data)
    min_sim = 0
    c_total = np.sum(total_data, axis = 0)
    for i, pixel in enumerate(total_data):
        c_total_i = c_total - total_data[i]
        #data_sets = [np.append(i_sum, len(total_data) - 1)]
        Indices = gen_sim_dict(c_total_i, n_fingerprints - 1)
        sim_index = Indices[weight][n_ary]
        if sim_index > min_sim:
            min_sim = sim_index
            index = i
        else:
            pass
    return index

#for m in [100000, 1000000, 10000000]:
#    c_total = np.array(list(range(m)), dtype='float32')
#    start = time.time()
#    counters = calculate_counters(c_total, n_fingerprints=2)
#    t = time.time() - start
#    print(t)

#total_data = np.array([[0,1,0,0,1,0],[1,0,1,1,0,1],[1,0,0,0,1,1],[1,1,0,1,1,1],[0,1,1,0,1,1]])

#fp_total = 10000
#fp_size = 100000
#total_data = np.random.randint(2, size=(fp_total, fp_size), dtype='int8')
#total_sum = np.sum(total_data, axis=0)
#
#old = 0
#
#if old:
#    start = time.time()
#    total = calculate_comp_sim(total_data, total_sum, c_threshold=None, n_ary = 'RR', weight = 'nw')
#    t = time.time() - start
#    print(t)
#else:
#    start = time.time()
#    total = calculate_comp_sim_rr(total_data, total_sum)
#    t = time.time() - start
#    print(t)
#    print(total)

#dim = 100
#total_data = np.random.random((dim, dim))
#c_total = real_var2_pre(total_data)
#c_threshold = 0.56
#d = gen_sim_dict(c_total, dim, c_threshold=c_threshold)
#print(d['nw']['RR'])
