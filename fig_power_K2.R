##------------------------------------------------------------------------------
## Imports libraries 

library(fastcluster)
library(clusterpval)
library(ggplot2)
library(patchwork)
library(latex2exp)
source("functions.R")

##------------------------------------------------------------------------------
## Specifies settings 

# settings for data generating process
ss = 1
n = 30
q = 2
vec_delta = c(0, 1, 2, 3, 4, 5, 6, 7)
ld = length(vec_delta)
cl_true = rep(1:2, each = round(n / 2))

# settings for clustering 
K = 2

# other settings 
thresh = 0.05 # significance level
num_trial = 500 # number of trials

##------------------------------------------------------------------------------
## Creates a data set consistent with the alternative for low signal strength

delta = vec_delta[4]

set.seed(1)
X = fun_gen_X(n, q, ss, delta)

# runs clustering
hc = hclust(dist(X) ** 2, method = "average")
cl = stats::cutree(hc, K)

# visualizes
df_low = data.frame(x = X[, 1],
                    y = X[, 2],
                    col_vec = as.factor(cl),
                    shape_vec = as.factor(rep(c("A", "B"),
                                              each = n / 2)))
plot_sim_low = ggplot(df_low, aes(x = x, y = y)) +
  geom_point(aes(shape = shape_vec, color = col_vec), size = 5, stroke = 3) +
  scale_color_manual(values = c("1" = "#898989",
                                "2" = "#64379E")) +
  scale_shape_manual(values = c("A" = 15,
                                "B" = 16)) +
  labs(title = TeX(paste0("Visualization for ", "$\\delta=$",
                          delta)),
       color = "Clusters Formed", shape = "True Clusters",
       x = NULL, y = "Setting 1") +
  theme(legend.position = "none") + 
  xlim(-2, 7) + 
  ylim(-3, 2) + 
  plot_details + 
  theme(axis.title = element_text(size = 40))

##------------------------------------------------------------------------------
## Creates a data set consistent with the alternative for high signal strength

delta = vec_delta[6]

set.seed(1)
X = fun_gen_X(n, q, ss, delta)

# runs clustering
hc = hclust(dist(X) ** 2, method = "average")
cl = stats::cutree(hc, K)

# visualizes
df_high = data.frame(x = X[, 1],
                     y = X[, 2],
                     col_vec = as.factor(cl),
                     shape_vec = as.factor(rep(c("A", "B"),
                                               each = n / 2)))
plot_sim_high = ggplot(df_high, aes(x = x, y = y)) +
  geom_point(aes(shape = shape_vec, color = col_vec), size = 5, stroke = 3) +
  scale_color_manual(values = c("1" = "#898989",
                                "2" = "#64379E")) +
  scale_shape_manual(values = c("A" = 15,
                                "B" = 16)) +
  labs(title = TeX(paste0("Visualization for ", "$\\delta=$",
                          delta)),
       color = "Clusters Formed", shape = "True Clusters",
       x = NULL, y = NULL) +
  theme(legend.position = "none") + 
  xlim(-2, 7) + 
  ylim(-3, 2) + 
  plot_details

##------------------------------------------------------------------------------
## Computes empirical power

mat_power = matrix(0, 4, ld)
mat_power_eb = matrix(0, 4, ld)

for (l in 1:4) {
  for (i in 1:ld) {
    delta = vec_delta[i]
    vec_pval = rep(0, num_trial)
    set.seed(1)
    for(j in 1:num_trial) {
      # generates data
      X = fun_gen_X(n, q, ss, delta)
      
      # runs clustering
      hc = hclust(dist(X) ** 2, method = "average")
      cl = stats::cutree(hc, K) 
      
      # checks if the pair k1, k2 is under the alternative
      ind_pair = c(which(cl == 1), which(cl == 2))
      if (length(unique(cl_true[ind_pair])) == 1) {
        vec_pval[j] = NA
      } else {
        # computes p-values
        if (l == 1) { 
          vec_pval[j] = fun_proposed_exact(X)
          
        } else if (l == 2) { 
          ss_hat_all = fun_ss_hat_all(X, n, q)
          vec_pval[j] = test_hier_clusters_exact(X, "average", hc, K, 1, 2, 
                                                 sig = sqrt(ss_hat_all))$pval
          
        } else if (l == 3) { 
          ss_hat_clustered = fun_ss_hat_clustered(X, n, q, cl, K)
          vec_pval[j] = test_hier_clusters_exact(X, "average", hc, K, 1, 2, 
                                                 sig = sqrt(ss_hat_clustered))$pval
          
        } else { 
          vec_pval[j] = test_hier_clusters_exact(X, "average", hc, K, 1, 2, 
                                                 sig = sqrt(ss))$pval
        }
      }
    }
    # computes power
    tmp = fun_power(vec_pval, thresh)
    mat_power[l, i] = tmp[1]
    mat_power_eb[l, i] = tmp[2]
  }
}

##------------------------------------------------------------------------------
## Creates plots of empirical power

stack = c(mat_power[1, ], mat_power[2, ], mat_power[3, ], mat_power[4, ])
stack_eb = c(mat_power_eb[1, ], mat_power_eb[2, ], mat_power_eb[3, ], 
             mat_power_eb[4, ])

df = data.frame(delta = rep(vec_delta, 4),  
                power = stack, 
                case = factor(rep(c("proposed", "all", 
                                    "clustered", "oracle"), each = ld)), 
                leb = stack - stack_eb, 
                ueb = stack + stack_eb)

plot_power = ggplot(df, aes(delta, power, group = case)) + 
  geom_errorbar(aes(ymin = leb, ymax = ueb, col = case), 
                width = 0.3, size = 2) + 
  geom_line(aes(col = case), linetype = "solid", size = 2) + 
  geom_point(aes(col = case), size = 2) + 
  labs(title = "Empirical Power", 
       x = TeX(r'($\delta$)'), y = "Empirical Power", color = NULL, 
       linetype = NULL) +
  theme(legend.position = "none", legend.justification = c(1, 0)) + 
  scale_color_manual(values = c("proposed" = "#C8375C",
                                "all" = "#3E9FB3", 
                                "clustered" = "#FFC918", 
                                "oracle" = "#605D58")) + 
  ylim(0, 1) + 
  plot_details

##------------------------------------------------------------------------------
## Saves workspace 

to_save = "workspaces/fig_power_K2.RData"
if (dir.exists("workspaces") == 0) {
  dir.create("workspaces")
  save.image(file = to_save)
} else {
  save.image(file = to_save)
}

##------------------------------------------------------------------------------
## Saves plot

plot_comb = plot_sim_low | plot_sim_high | plot_power
to_save = paste0("plots/fig_power_K2.png")
if (dir.exists("plots") == 0) {
  dir.create("plots")
  ggsave(plot_comb, file = to_save, width = 26, height = 6.5)
} else {
  ggsave(plot_comb, file = to_save, width = 26, height = 6.5)
}

##------------------------------------------------------------------------------