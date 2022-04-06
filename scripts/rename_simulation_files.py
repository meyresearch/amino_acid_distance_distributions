"""Code to rename files from old repo."""
import glob
import os

path_to_sim_files = f"../../data/simulation_data_matrices/"
sim_file_lengths = ["100", "209", "302"]
correct_lengths = ["100", "200", "300"]

for i in range(len(sim_file_lengths)):
    files = glob.glob(path_to_sim_files+f"Nora50_matrix_{sim_file_lengths[i]}_*")
    for j in range(len(files)):
        bash = f"cp {files[j]} ../data/simulations/3d/matrices/matrix_{correct_lengths[i]}_{j}.txt"
        os.system(bash)
