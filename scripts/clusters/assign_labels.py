import numpy as np
import sys
sys.path.insert(0, "../../")
from sklearn.cluster import KMeans
import os

# System info - EDIT THESE
input_traj_numpy = '../../example/aligned_tau.npy'
N_atoms = 50
sieve = 1

# K-means params - EDIT THESE
n_clusters = 20
output_dir = 'outputs'

if __name__ == '__main__':
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    traj_numpy = np.load(input_traj_numpy)[::sieve]
    mod = KMeans(n_clusters=n_clusters, n_init=10, random_state=0)
    labels, centers, n_iter = mod.fit(traj_numpy).labels_, mod.cluster_centers_, mod.n_iter_
    sort_labels_by_size = np.argsort(np.bincount(labels))[::-1]
    labels = np.array([np.where(sort_labels_by_size == i)[0][0] for i in labels])
    
    # Save cluster labels
    with open(f'{output_dir}/labels_{n_clusters}.csv', 'w') as f:
        f.write('# Frame Index, Cluster Index\n')
        for i, row in enumerate(labels):
            f.write(f'{i * sieve},{row}\n')
    
    # Calculate population of each cluster
    with open(f'{output_dir}/summary_{n_clusters}.txt', 'w') as f:
        f.write(f'# Cluster Index, Fraction out of total pixels {len(labels)}\n')
        for i, row in enumerate(np.bincount(labels)):
            f.write(f'{i} {row}\n')