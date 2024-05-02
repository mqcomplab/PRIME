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
    - [2. NANI Screening](#2-nani-screening)
    - [3. Analysis of NANI Screening Results](#3-analysis-of-nani-screening-results)
    - [4. Cluster Assignment](#4-cluster-assignment)
    - [5. Extract frames for each cluster (Optional)](#5-extract-frames-for-each-cluster-optional)

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
```
`modules` contains the functions required to run the algorithm. `sample_clusters/clusttraj.c*` contains sample clustering files prepared through CPPTRAJ Hierarchical clustering. `new_clusters/normalize.py` normalizes clustering files through min-max normalization. `similarity.py` generates a similarity dictionary from running the protein refinement method.

## Tutorial

### Step 1.
```
git clone https://github.com/lexin-chen/PRIME.git
```
### Step 2. Normalization

- Normalize the trajectory data between $[0,1]$ using the Min-Max Normalization. 
```
cd sample_clusters
python normalize.py -b 950
```
Output will have `normed_clusttraj.c*` files and `normed_data.txt`, appended all normed files together.

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
