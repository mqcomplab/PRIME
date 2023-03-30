import sys
sys.path.insert(0, '../../modules')
from sim_modules_vector import *
from modules import *
import numpy as np

break = 950
# Combine all the files into one files,
# Normalize all data using the min-max normalization,
# and getting the min and max for the whole dataset.

data = read_cpptraj(break)
norm = Normalize(data=data)
normed_data = norm.get_normed_data()
min, max = norm.get_min_max()
np.savetxt("normed_data.txt", normed_data)

# Normalize each file using the min and max determined from above. 
data = read_cpptraj(break, min=min, max=max, normalize=True)
