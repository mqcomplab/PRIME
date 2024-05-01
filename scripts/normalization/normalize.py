import numpy as np
import sys
sys.path.insert(0, '../../../../../clustering/alignedClustering/')
from src.inputs.preprocess import gen_traj_numpy, normalize_file, Normalizer
import re
import glob

# System info - EDIT THESE
input_top = '../../topol.top'
output_base_name = 'normed_clusttraj'
atomSelection = 'resid 2:17 21:28 34:38 42:58 64:69 73:77 and name CA CB'

if __name__ == '__main__':
    list_clusttraj = sorted(glob.glob("../outputs/clusttraj_*"), key=lambda x: int(re.findall("\d+", x)[0]))
    list_clusttraj = list_clusttraj[:11]
    all_clusttraj = []
    for clusttraj in list_clusttraj:
        traj_numpy = gen_traj_numpy(input_top, clusttraj, atomSelection)
        all_clusttraj.append(traj_numpy)
        print(traj_numpy.shape)
    print(len(all_clusttraj))
    concat_clusttraj = np.concatenate(all_clusttraj)
    normed_data, min, max, avg = normalize_file(concat_clusttraj, norm_type='v3')
    np.savetxt('normed_data.txt', normed_data)
    print(f'Min: {min}, Max: {max}, Avg: {avg}')
    for i, traj in enumerate(all_clusttraj):
        norm = Normalizer(data=traj, custom_min=min, custom_max=max)
        normed_frame = norm.get_v3_norm()
        np.savetxt(f'{output_base_name}.c{i}', normed_frame)
    