# Protein Refinement in Molecular Ensembles (PRIME)
**🪄 Mastering the art of protein refinement, one conformation at a time 🪄**
Refining protein structures is important because the accuracy of protein structures influence how our understanding of its function and its interactions with other molecules, which can help to design new drugs to target specific molecular interactions.  

This repo contains six different ways of determining the native structure of biomolecules from simulation or clustering data. 

### Usage
`modules` contains the functions required to run the algorithm. `scripts` contains many examples of how to run the protein refinery framework.

Step 1. 
- Normalize the trajectory data between $[0,1]$ using the Min-Max Normalization. 
- In `sample_clusters/spacing.py`, edit `break =` for the number of lines of one frame. 
- This will give `normed_data.txt` and `normed_*.txt` for each file. 

Step 2. Execute the scripts in `scripts`
- 2.1 refers to the pairwise method.
- 2.2 refers to the union method.
- 3.1 refers to the medoid method.
- 3.2 refers to the outlier method. 

