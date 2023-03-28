import sys
sys.path.insert(0, '../../modules')
from sim_modules_vector import *
from modules import *
import numpy as np

data = read_cpptraj(950)
norm = Normalize(data=data)
normed_data = norm.get_normed_data()
min, max = norm.get_min_max()
np.savetxt("normed_data.txt", normed_data)
#with open("minmax.txt", "w") as f:
#	f.write(f"min is {min}\n")
#	f.write(f"max is {max}")

data = read_cpptraj(950, min=min, max=max, normalize=True)

