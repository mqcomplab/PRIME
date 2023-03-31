<img src="img/header.jpg" width="1000" height=auto align="center"></a>

**🪄 Mastering the art of protein refinement, one conformation at a time 🪄**

Refining protein structures is important because the accuracy of protein structures influence how our understanding of its function and its interactions with other molecules, which can help to design new drugs to target specific molecular interactions.  

This repo contains six different ways of determining the native structure of biomolecules from simulation or clustering data. 

<img src="img/methods.jpg" alt="Girl in a jacket" width="500" height=auto align="center"></a>

*Fig 1. Six techniques of protein refinement. Blue is top cluster.* 

### Usage
`modules` contains the functions required to run the algorithm. `scripts` contains many examples of how to run the protein refinery framework. The naming schemes followes the *Fig. 1*. 

Step 1. 
- Normalize the trajectory data between $[0,1]$ using the Min-Max Normalization. 
- In `sample_clusters/spacing.py`, edit `break =` for the number of lines of one frame. 
- This will give `normed_data.txt` and `normed_*.txt` for each file. 



Step 2. Execute the scripts in `scripts`. Naming schemes is `{method}_{indice}`. The default indice is Russell Rao. For example `2.1.py` uses 2.1 method from Figure on the right using Russell Rao index. `2.1_sm.py` uses 2.1 method using the Sockal Michener index. All the scripts will generate a dictionary of similarity results.
The dictionary format is 
```
"""
Keyes are frame #.
Values are [cluster 1 similarity, cluster #2 similarity, 
            ..., average similarity of all clusters]
"""
{f1: [0.1,
      0.2,
      ...
      average],
 f2: [0.2,
      0.3,
      ...
      average],
```
Step 3. Determine the native structure using `scripts/rep.py`.

