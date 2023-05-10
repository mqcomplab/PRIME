import modules as mod

mod.graph_rep_frames(folder_pattern="[Nn]*", sim_folder="nw", trim_frac=None)
mod.graph_rep_frames(folder_pattern="[Nn]*", sim_folder="nw", trim_frac=0.1)
mod.graph_rep_frames(folder_pattern="[Nn]*", sim_folder="nw", trim_frac=0.2)
mod.graph_rep_frames(folder_pattern="[Nn]*", sim_folder="nw", trim_frac=0.3)

a = mod.graph_rep_frames(folder_pattern="[Nn]*", sim_folder="nw", n_ary="SM", trim_frac=None, return_dict=True)
b = mod.graph_rep_frames(folder_pattern="[Nn]*", sim_folder="nw", n_ary="SM", trim_frac=0.1, return_dict=True)
c = mod.graph_rep_frames(folder_pattern="[Nn]*", sim_folder="nw", n_ary="SM", trim_frac=0.2, return_dict=True)
d = mod.graph_rep_frames(folder_pattern="[Nn]*", sim_folder="nw", n_ary="SM", trim_frac=0.3, return_dict=True)
