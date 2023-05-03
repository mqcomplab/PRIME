import numpy as np
from modules.sim_modules_vector import *
from itertools import chain
import re
import glob

class Normalize:
    def __init__(self, file_path=None, data=None, custom_min=None, custom_max=None):
        if file_path is not None:
            self.file_path = file_path
            self.data = np.genfromtxt(self.file_path)
        elif data is not None:
            self.data = data
        if custom_min is not None and custom_max is not None:
            self.min = custom_min
            self.max = custom_max
        else:
            self.min = np.min(self.data)
            self.max = np.max(self.data)
        self.normed_data = (self.data - self.min) / (self.max - self.min)
        self.esim_norm = 1 - np.abs(self.normed_data - np.mean(self.normed_data, axis=0))
        #self.c_total = np.sum(1 - np.abs(self.normed_data - np.mean(self.normed_data, axis=0)), axis=0)
        self.c_total = np.sum(self.normed_data, axis=0)
    def get_c_total(self):
        return self.c_total
    def get_normed_data(self):
        return self.normed_data
    def get_min_max(self):
        return self.min, self.max
    def get_esim_norm(self):
        return self.esim_norm

def read_cpptraj(break_line, min=None, max=None, normalize=False):
    """ Read CPPTRAJ files to ensure proper spacings for future post-processing """
    input_files = sorted(glob.glob("clusttraj.c*"), key=lambda x: int(re.findall("\d+", x)[0]))
    break_line = break_line
    frames_list = []
    count_frames = []
    for file in input_files:
        with open(file, 'r') as infile:
            lines = [line.rstrip() for line in infile][1:]
        sep_lines = [[line[i:i+8] for i in range(0, len(line), 8)] for line in lines]
        chunks = [sep_lines[i:i+break_line] for i in range(0, len(sep_lines), break_line)]
        str_frames = [list(chain.from_iterable(chunk)) for chunk in chunks]
        str_frames = [' '.join(frame) for frame in str_frames]
        frames = np.array([np.fromstring(frame, dtype='float32', sep=' ') for frame in str_frames])
        if normalize is True:
            norm = Normalize(data=frames, custom_min=min, custom_max=max)
            normed_frame = norm.get_normed_data()
            np.savetxt(f"normed_{file}", normed_frame)
        else:
            frames_list.append(frames)
        count_frames.append(len(frames))
    if normalize is False:
        data = np.concatenate(frames_list, axis=0)
        return data
    