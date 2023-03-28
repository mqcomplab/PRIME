def write_cpptraj_script(cluster_folder, atom_sel):
    # Define folder path and atom selection using CPPTRAJ syntax.
    topol_file = "topol.top"
    traj_files = [f"{cluster_folder}/clusttraj.c{i}" for i in range(0, 11)]
    ref_traj = f"{cluster_folder}/clusttraj.c0"

    # Read frame number from file
    with open("centroid_frame_number.dat", "r") as f:
        line = f.readline().strip()
        frame_num = int(line.split()[-1])

    # Build cpptraj command
    cpptraj_command = f"parm {topol_file}\n"
    cpptraj_command += "\n".join([f"trajin {traj_file}" for traj_file in traj_files])
    cpptraj_command += f"\nreference {ref_traj} {frame_num}"
    cpptraj_command += f"\nrmsd protein_rmsd {atom_sel} reference out rmsd.dat\n"
    cpptraj_command += "run"

    # Write cpptraj command to file
    with open("rmsd.in", "w") as f:
        f.write(cpptraj_command)

