<img src="img/3.png" width="1000" height=auto align="center"></a>
<br>

<h3 align="center">
    <p><b>🪄 Mastering the art of protein refinement, one conformation at a time 🪄</b></p>
    </h3>

Refining protein structures is important because the accuracy of protein structures influence how our understanding of its function and its interactions with other molecules, which can help to design new drugs to target specific molecular interactions. This repo contains six different ways of determining the native structure of biomolecules from simulation or clustering data. 

<img src="img/methods.jpg" alt="Girl in a jacket" width="500" height=auto align="center"></a>

*Fig 1. Six techniques of protein refinement. Blue is top cluster.* 

## Usage
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
- `-m` - methods (*required)
- `-n` - number of clusters (*required)
- `-i` - similarity index (*required)
- `-t` - Fraction of outliers to trim in decimals. 
- `-w` - Weighing clusters by frames it contains.
- `-d` - directory where the `normed_clusttraj.c*` files are located.
- `-s` - location where `summary` file is located

```
python similarity.py -m medoid -n 11 -i RR
```
To generate a similarity dictionary using data in `sample_clusters/normed_clusttraj.c*` using the medoid method (3.1 in *Fig 1*) and Russell Rao index.

The result is a dictionary organized as followes:
Keys are frame #. Values are [cluster 1 similarity, cluster #2 similarity, ..., average similarity of all clusters].

### Step 4. Determine the native structure.

