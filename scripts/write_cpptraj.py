import sys
sys.path.insert(0, '../')
import modules as mod

cluster_folder = "Cpptraj_linkage_sieve_eps_all_ss_1"
atom_sel = ":3,4,5,6,7,14,15,16,17,18,19,29,30,31,40,41,42,43,45,46,51,52,53,54,55,56,57,60,61,62,63,64,65,71,72,73,75,76,77,78,79,80,81@CA,CB"
mod.write_cpptraj_script(cluster_folder, atom_sel)
