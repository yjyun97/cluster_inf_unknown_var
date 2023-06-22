This code reproduces the experiments in "Selective inference for clustering with unknown variance" (authors: Youngjoo Yun and Rina Foygel Barber), available at [http://arxiv.org/abs/2301.12999](http://arxiv.org/abs/2301.12999). 

Each of the files starting with "fig_" corresponds to a figure in the paper, and the correspondence is shown below. Running each of these files saves the same figure into a folder called "plots" and saves the workspace into a folder called "workspaces" in the working directory. Running **fig_penguins.R** further saves the results in [Table 1] into a folder called "results". 

* [Figure 1] **fig_motivation.R**
* [Figure 2] **fig_decomposition.R**
* [Figure 3] **fig_type1.R**
* [Figure 4] **fig_power_K2.R**, **fig_power_K3_equidistant.R**, and **fig_power_K3_horizontal.R**
* [Figure 5] **fig_type1_t5.R**
* [Figure 6] **fig_type1_t10.R**
* [Figure 7] **fig_type1_noniso.R**
* [Figure 8] **fig_penguins.R**

The functions used in the scripts above are stored in the file called **functions.R**.