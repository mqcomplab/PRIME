<img src="img/logo.png" width="800" height=auto align="center"></a>
<br>
## Protein Retrieval via Integrative Molecular Ensembles (PRIME)

<h3 align="center"> 
    <p><b>🪄 Predict Protein Structure with Precision 🪄</b></p>
    </h3>
Table of Contents
=================
- [Overview](#overview)
- [Tutorial](#tutorial)
    - [1. Input Preparations](#1-input-preparations)
    - [2. Cluster Assignment](#2-cluster-assignment)
    - [3. Cluster Normalization](#3-cluster-normalization)

## Overview
<p>Protein structures prediction is important because the accuracy of protein structures influence how our understanding of its function and its interactions with other molecules, which can help to design new drugs to target specific molecular interactions. This repo contains six different ways of determining the native structure of biomolecules from simulation or clustering data. <b>These methods perfectly mapped all the structural motifs in the studied systems and required unprecedented linear scaling.</b></p>
&nbsp
<figure>
<img src="img/2k2e.png" alt="2k2e" width="400" height=auto align="right"></a>
    <figcaption><i>Fig 1. Superposition of the most representative structures found with extended indices (yellow) and experimental native structures (blue) of 2k2e.</i></figcaption>
</figure>

## Installation
``` 
git clone https://github.com/lexin-chen/PRIME.git
cd PRIME
```
requires Python 3.6+ and the following packages: MDAnalysis, numpy, and matplotlib. 

`modules` contains the functions required to run the algorithm. `sample_clusters/clusttraj.c*` contains sample clustering files prepared through CPPTRAJ Hierarchical clustering. `new_clusters/normalize.py` normalizes clustering files through min-max normalization. `similarity.py` generates a similarity dictionary from running the protein refinement method.

## Tutorial
The following tutorial will guide you through the process of determining the native structure of a biomolecule using the PRIME algorithm. If you already have clustered data, you can skip to Step 3.

### 1. Input Preparations
Preparation for Molecular Dynamics Trajectory

Prepare a valid topology file (e.g. `.pdb`, `.prmtop`), trajectory file (e.g. `.dcd`, `.nc`), and the atom selection. This step will convert a Molecular Dynamics trajectory to a numpy ndarray. **Make sure the trajectory is already aligned and/or centered if needed!**

**Step-by-step tutorial can be found in the [scripts/inputs/preprocessing.ipynb](../scripts/inputs/preprocessing.ipynb).**

### 2. Cluster Assignment
In this example, we will use *k*-means clustering to assign labels to the clusters and the number of clusters will be 20. Any clustering method can be used as long as the data is clustered (e.g. DBSCAN, Hierarchical Clustering). **Please check out [MDANCE](https://github.com/mqcomplab/MDANCE) for more clustering methods!**

[scripts/nani/assign_labels.py](scripts/nani/assign_labels.py) will assign labels to the clusters using *k*-means clustering

    # System info - EDIT THESE
    input_traj_numpy = '../../example/aligned_tau.npy'
    N_atoms = 50
    sieve = 1

    # K-means params - EDIT THESE
    n_clusters = 20
    output_dir = 'outputs'

#### Inputs
##### System info
`input_traj_numpy` is the numpy array prepared from step 1. <br>
`N_atoms` is the number of atoms used in the clustering, should be same as atom selection in step 1.<br>
`sieve` takes every `sieve`th frame from the trajectory for analysis. <br>

##### *k*-means params
`n_clusters` is the number of clusters for labeling. <br>
`output_dir` is the directory where the output files will be saved. <br>

#### Execution
```bash
python assign_labels.py
```
#### Outputs
1. csv file containing the cluster labels for each frame.
2. csv file containing the population of each cluster.

### 3. Cluster Normalization
With already clustered data, [scripts/normalization/normalize.py](scripts/normalization/normalize.py) Normalize the trajectory data between $[0,1]$ using the Min-Max Normalization. 

    # System info - EDIT THESE
    input_top = '../../example/aligned_tau.pdb'
    unnormed_cluster_dir = '../clusters/outputs/clusttraj_*'
    output_base_name = 'normed_clusttraj'
    atomSelection = 'resid 3 to 12 and name N CA C O H'
    n_clusters = 10

#### Inputs
##### System info
`input_top` is the topology file used in the clustering. <br>
`unnormed_cluster_dir` is the directory where the clustering files are located from step 2. <br>
`output_base_name` is the base name for the output files. <br>
`atomSelection` is the atom selection used in the clustering. <br>
`n_clusters` is the number of clusters used in the PRIME. If number less than total number of cluster, it will take top *n* number of clusters. <br>

```bash
python normalize.py
```

#### Outputs
1. `normed_clusttraj.c*.npy` files, normalized clustering files.
2. `normed_data.npy`, appended all normed files together.

### Step 3. Similarity Calculations
Back at root directory, `similarity.py` generates a similarity dictionary from running the protein refinement method. 
- `-h` - for help with the argument options.
- `-m` - methods (*required*)
- `-n` - number of clusters (*required*)
- `-i` - similarity index (*required*)
- `-t` - Fraction of outliers to trim in decimals. 
- `-w` - Weighing clusters by frames it contains.
- `-d` - directory where the `normed_clusttraj.c*` files are located.
- `-s` - location where `summary` file is located

```
python similarity.py -m medoid -n 11 -i RR -t 0.1
```
To generate a similarity dictionary using data in `sample_clusters/normed_clusttraj.c*` using the medoid method (3.1 in *Fig 1*) and Russell Rao index. In addition, 10% of the outliers were trimmed. 

The result is a dictionary organized as followes:
Keys are frame #. Values are [cluster 1 similarity, cluster #2 similarity, ..., average similarity of all clusters].

### Step 4. Determine the native structure.
`scripts/rep.py`

### Further Reading
<img src="img/methods.jpg" alt="methods" width="500" height=auto align="center"></a>

 *Fig 2. Six techniques of protein refinement. Blue is top cluster.* 
