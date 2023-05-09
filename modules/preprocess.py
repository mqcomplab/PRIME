import numpy as np
from itertools import chain
import re
import glob

class Normalize:
    """Class for normalizing data.

    Args:
    file_path (str, optional): The file path of the input data. If provided, the data is read from
                               the file. Defaults to None.
    data (np.ndarray, optional): The input data as a numpy array. If provided, the file_path argument
                                  is ignored. Defaults to None.
    custom_min (float or None, optional): The minimum value to use for normalization. If not provided, 
                                          the minimum value of the input data is used. Defaults to None.
    custom_max (float or None, optional): The maximum value to use for normalization. If not provided, 
                                          the maximum value of the input data is used. Defaults to None.

    Attributes:
    normed_data (np.ndarray): The normalized input data as a numpy array.
    c_total (np.ndarray): The sum of the absolute difference of each data point from the mean along the 
                          column axis of the normalized data.
    esim_norm (np.ndarray): The element-wise similarity between each data point and the mean along the 
                            column axis of the normalized data.
    min (float): The minimum value of the input data.
    max (float): The maximum value of the input data.
    """    
    def __init__(self, file_path=None, data=None, custom_min=None, custom_max=None):
        if file_path:
            self.file_path = file_path
            self.data = np.genfromtxt(self.file_path)
        elif data is not None:
            self.data = data
        if custom_min and custom_max:
            self.min = custom_min
            self.max = custom_max
        else:
            self.min = np.min(self.data)
            self.max = np.max(self.data)
        self.v3_norm = (self.data - self.min) / (self.max - self.min)
        self.v2_norm = 1 - np.abs(self.v3_norm - np.mean(self.v3_norm, axis=0))
        self.c_total = np.sum(1 - np.abs(self.v3_norm - np.mean(self.v3_norm, axis=0)), axis=0)
    
    def get_min_max(self):
        return self.min, self.max
    
    def get_v2_norm(self):
        return self.v2_norm
    
    def get_v3_norm(self):
        return self.v3_norm
    
    def get_c_total(self):
        return self.c_total

def read_cpptraj(break_line=None, norm_type=None, min=None, max=None, normalize=False):
    """Read CPPTRAJ files to convert to numpy ndarray formatting and normalize the data.
    
    Args:
    break_line (int): The number of columns per line of the input file. (have to n-1 because ignore first line)
    min (float or None, optional): The minimum value to use for normalization. If not provided, 
                                    the minimum value of the input data is used. Defaults to None.
    max (float or None, optional): The maximum value to use for normalization. If not provided, 
                                    the maximum value of the input data is used. Defaults to None.
    normalize (bool, optional): Whether to normalize the input data. If True, the data is 
                                normalized to the range [0, 1]. Defaults to False.
                              
    Returns:
    np.ndarray: The concatenated input data as a numpy array.
    """
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
        if normalize:
            norm = Normalize(data=frames, custom_min=min, custom_max=max)
            if norm_type == "v2":
                normed_frame = norm.get_v2_norm()
            elif norm_type == "v3":
                normed_frame = norm.get_v3_norm()
            np.savetxt(f"normed_{file}", normed_frame)
        else:
            frames_list.append(frames)
        count_frames.append(len(frames))
    if not normalize:
        data = np.concatenate(frames_list, axis=0)
        return data