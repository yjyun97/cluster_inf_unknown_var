This code reproduces the experiments in "Selective inference for clustering with unknown variance" (Youngjoo Yun and Rina Foygel Barber), available at [http://arxiv.org/abs/2301.12999](http://arxiv.org/abs/2301.12999). 

Each of the files starting with "sim_" (which stands for simulations) and "real_" (which stands for real data application) corresponds to a figure in the paper, and the correspondence is shown below. Running each of these files saves the same figure into a folder called "plots" and saves the workspace into a folder called "workspaces" in the working directory. Running real_penguins.R  further saves the results into a folder called "results". 

[Figure 1] sim_motivation.R
[Figure 2] sim_decomposition.R
[Figure 3] sim_type1.R
[Figure 4] sim_power_K_2.R, sim_power_K_3_equidistant.R, and sim_power_K_3_horizontal.R
[Figure 5] real_penguins.R

The functions used in the scripts above are stored in the file called functions.R.
