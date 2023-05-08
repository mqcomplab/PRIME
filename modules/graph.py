import numpy as np
import matplotlib.pyplot as plt
import glob
import re
import matplotlib
matplotlib.use('Agg')
import os

def graph_rep_frames(folder_pattern="[Nn]*",weighted=True, n_ary="RR", trim_frac=0.1, return_dict=None):
    if not os.path.exists("graphs"):
        os.makedirs("graphs")
    input_folders = sorted(glob.glob(folder_pattern), key=lambda x: int(re.findall("\d+", x)[0]))
    if weighted is True:
        w = "w_"
    elif weighted is False:
        w = ""
    if trim_frac is not None:
        t = f"_t{int(float(trim_frac) * 100)}"
    elif trim_frac is None:
        t= ""
    
    # Dictionary for all RMSD results for each method.
    data_dict = {"medoid_all": [], "medoid_c0": [], "medoid_c0_trimmed": [], "pairwise": [], "union": [], "medoid": [], "outlier": []}
    for folder in input_folders:
        data = np.loadtxt(f"{folder}/rmsd.dat", usecols=1, skiprows=1)
        list = np.genfromtxt(f"{folder}/w/{w}rep_{n_ary}{t}.txt", delimiter=",")
        selected_lines = list.astype(int).tolist()
        for i, key in enumerate(data_dict.keys()):
            data_dict[key].append(data[selected_lines[i]])
    
    # Create new dictionary for rmsd(method) - rmsd(med_c0) values.
    n_dict = {"medoid_all": [], "medoid_c0_trimmed": [], "pairwise": [], "union": [], "medoid": [], "outlier": []}
    for key in ["medoid_all", "medoid_c0_trimmed", "pairwise", "union", "medoid", "outlier"]:
        n_dict[key] = [a_i - b_i for a_i, b_i in zip(data_dict[key], data_dict["medoid_c0"])]
    d_means = []
    for k in n_dict:
        average = sum(n_dict[k]) / len(n_dict[k])
        d_means.append(average)
    
    # Plotting points for each system
    x_values = [r"$med_{c_{all}}$", r"$med_{c_{0}~trimmed}$", r"$\langle sim(F_{c0}, F_k)\rangle$", r"$esim({F_{c0}}\cup C_k)$", r"$sim(F_{C_0}, med_{C_k})$", r"$sim(F_{C_0}, out_{C_k})$"]
    plt.ylabel(r"$\mathrm{\Delta \langle RMSD\rangle~(Å)}$")
    new_ys = [n_dict["medoid_all"], n_dict["medoid_c0_trimmed"], n_dict["pairwise"], n_dict["union"], n_dict["medoid"], n_dict["outlier"]]
    for x, y in zip(x_values, new_ys):
        for yi in y:
            c = '#00b5c5' if yi < 0 else '#c57300'
            plt.scatter(x, yi, c=c)
    plt.scatter(x_values, d_means, marker="o", color="#0047AB", label="average", s=70)
    plt.axhline(y = 0, linestyle="--", color="black", zorder=0)
    plt.ylim(-5, 15)
    plt.legend()
    plt.savefig(f"graphs/{w}protein_{n_ary}{t}.png", bbox_inches="tight", dpi=300, transparent=True)
    plt.close()
    
    # Returning data dictionary if requested.
    if return_dict is not None:
        return n_dict

def graph_rep_frames_fracs(folder_pattern="[Nn]*",weighted=True, n_ary="RR", trim_frac=0.1, return_dict=None):
    if not os.path.exists("graphs"):
        os.makedirs("graphs")
    input_folders = sorted(glob.glob(folder_pattern), key=lambda x: int(re.findall("\d+", x)[0]))
    if weighted is True:
        w = "w_"
    elif weighted is False:
        w = ""
    if trim_frac is not None:
        t = f"_t{int(float(trim_frac) * 100)}"
    elif trim_frac is None:
        t= ""
    
    # Dictionary for all RMSD results for each method.
    data_dict = {"medoid_all": [], "medoid_c0": [], "pairwise": [], "union": [], "medoid": [], "outlier": []}
    for folder in input_folders:
        data = np.loadtxt(f"{folder}/rmsd.dat", usecols=1, skiprows=1)
        list = np.genfromtxt(f"{folder}/{w}rep_{n_ary}{t}.txt", delimiter=",")
        selected_lines = list.astype(int).tolist()
        for i, key in enumerate(data_dict.keys()):
            data_dict[key].append(data[selected_lines[i]])
    
    # Create new dictionary for rmsd(method) - rmsd(med_c0) values.
    n_dict = {"medoid_all": [], "pairwise": [], "union": [], "medoid": [], "outlier": []}
    for key in ["medoid_all", "pairwise", "union", "medoid", "outlier"]:
        n_dict[key] = [a_i - b_i for a_i, b_i in zip(data_dict[key], data_dict["medoid_c0"])]
    d_means = []
    fracs = []
    for k in n_dict:
        average = sum(n_dict[k]) / len(n_dict[k])
        d_means.append(average)
        num_points = len(n_dict[k])
        count_below_zero = 0
        for point in n_dict[k]:
            if point < 0:
                count_below_zero += 1
        fracs.append(count_below_zero / num_points)
    
    # Plot scatter plot for each system
    x_values = [r"$med_{c_{all}}$", r"$\langle sim(F_{c0}, F_k)\rangle$", r"$esim({F_{c0}}\cup C_k)$", r"$sim(F_{C_0}, med_{C_k})$", r"$sim(F_{C_0}, out_{C_k})$"]
    fig1, ax1 = plt.subplots()
    ax1.set_ylabel(r"$\mathrm{\Delta \langle RMSD\rangle~(Å)}$")
    new_ys = [n_dict["medoid_all"], n_dict["pairwise"], n_dict["union"], n_dict["medoid"], n_dict["outlier"]]
    for x, y in zip(x_values, new_ys):
        for yi in y:
            c = '#00ab64' if yi < 0 else '#c57300'
            ax1.scatter(x, yi, c=c)
    ax1.scatter(x_values, d_means, marker="o", color="#ab0047", s=70)
    ax1.axhline(y = 0, linestyle="--", color="black", zorder=0)
    ax1.set_ylim(-5, 15)

    labels = {'#ab0047': 'Average', '#00ab64': 'Below zero', '#c57300': 'Above zero'}
    handles = [plt.plot([], [], marker="o", ms=10, ls="", mec=None, color=color, label=label)[0] for color, label in labels.items()]
    ax1.legend(handles=handles, numpoints=1, loc='upper right', bbox_to_anchor=(1.0, 1.0), borderaxespad=0.1)

    # Plot bar graph of fraction of point crossed 0
    ax2 = ax1.twinx()
    ax2.set_ylabel("Fraction below zero")
    ax2.bar(x_values, fracs, width=0.4, color="#0047ab")
    #ax2.yaxis.set_major_locator(ticker.MultipleLocator(1 / len(n_dict["medoid"])))
    #ax2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, pos: str(Fraction(y).limit_denominator())))
    ax2.set_ylim(0, 1)
    ax1.set_zorder(ax2.get_zorder()+1)
    ax1.set_frame_on(False)
    fig1.savefig(f"graphs/{w}frac_{n_ary}{t}.png", bbox_inches="tight", dpi=300, transparent=True)

    # Returning data dictionary if requested.
    if return_dict is not None:
        return n_dict
