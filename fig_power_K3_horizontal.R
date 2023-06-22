##------------------------------------------------------------------------------
## Imports libraries 

library(fastcluster)
library(clusterpval)
library(ggplot2)
library(patchwork)
library(MASS)
library(latex2exp)
source("functions.R")

##------------------------------------------------------------------------------
## Specifies settings 

# settings for data generating process
ss = 1
ss_diag = ss * c(1, 1)
n = 30
n_each = round(n / 3)
q = 2
vec_delta = c(0, 1, 2, 3, 4, 5, 6, 7)
ld = length(vec_delta)
cl_true = rep(1:3, each = n_each)

# settings for clustering 
K = 3
k1 = 1
k2 = 2

# settings for importance sampling
ndraws = 8000 # number of samples drawn
alpha = 0.05 # sd of the proposal distribution

# other settings 
thresh = 0.05 # significance level
num_trial = 500 # number of trials

##------------------------------------------------------------------------------
## Creates a data set consistent with the alternative for low signal strength

delta = vec_delta[4]

set.seed(1)

mu1 = c(0, 0)
mu2 = c(delta, 0)
mu3 = c(2 * delta, 0)

X1 = mvrnorm(n = n_each, mu1, diag(ss_diag))
X2 = mvrnorm(n = n_each, mu2, diag(ss_diag))
X3 = mvrnorm(n = n_each, mu3, diag(ss_diag))

X = rbind(X1, X2, X3)

# runs clustering
hc = hclust(dist(X) ** 2, method = "average")
cl = stats::cutree(hc, K) 

# visualizes
df_low = data.frame(x = X[, 1],
                y = X[, 2],
                col_vec = as.factor(cl), 
                shape_vec = as.factor(rep(c("A", "B", "C"), 
                                          each = n_each)))

plot_sim_low = ggplot(df_low, aes(x = x, y = y)) + 
  geom_point(aes(shape = shape_vec, color = col_vec), size = 5, stroke = 3) + 
  scale_color_manual(values = c("1" = "#898989", 
                                "2" = "#64379E", 
                                "3" = "#67B3C9")) + 
  scale_shape_manual(values = c("A" = 15, 
                                "B" = 16, 
                                "C" = 17)) + 
  labs(title = TeX(paste0("Visualization for ", "$\\delta=$",
                          delta)), 
       color = "Clusters Formed", shape = "True Clusters", 
       x = NULL, y = "Setting 3") +
  theme(legend.position = "none") +
  xlim(-2, 12) + 
  ylim(-2, 2) + 
  plot_details + 
  theme(axis.title = element_text(size = 40))

##------------------------------------------------------------------------------
## Creates a data set consistent with the alternative for high signal strength

delta = vec_delta[6]

set.seed(1)

mu1 = c(0, 0)
mu2 = c(delta, 0)
mu3 = c(2 * delta, 0)

X1 = mvrnorm(n = n_each, mu1, diag(ss_diag))
X2 = mvrnorm(n = n_each, mu2, diag(ss_diag))
X3 = mvrnorm(n = n_each, mu3, diag(ss_diag))

X = rbind(X1, X2, X3)

# runs clustering
hc = hclust(dist(X) ** 2, method = "average")
cl = stats::cutree(hc, K) 

# visualizes
df_high = data.frame(x = X[, 1],
                y = X[, 2],
                col_vec = as.factor(cl), 
                shape_vec = as.factor(rep(c("A", "B", "C"), 
                                          each = n_each)))

plot_sim_high = ggplot(df_high, aes(x = x, y = y)) + 
  geom_point(aes(shape = shape_vec, color = col_vec), size = 5, stroke = 3) + 
  scale_color_manual(values = c("1" = "#898989", 
                                "2" = "#64379E", 
                                "3" = "#67B3C9")) + 
  scale_shape_manual(values = c("A" = 15, 
                                "B" = 16, 
                                "C" = 17)) + 
  labs(title = TeX(paste0("Visualization for ", "$\\delta=$",
                          delta)), 
       color = "Clusters Formed", shape = "True Clusters", 
       x = NULL, y = NULL) +
  xlim(-2, 12) + 
  ylim(-2, 2) + 
  theme(legend.position = "none") +
  plot_details

##------------------------------------------------------------------------------
## Computes empirical power 

mat_power = matrix(0, 4, ld) 
mat_power_eb = matrix(0, 4, ld)

for (i in 1:ld) {
  delta = vec_delta[i]
  set.seed(2)
  vec_pval1 = rep(0, num_trial)
  vec_pval2 = rep(0, num_trial)
  vec_pval3 = rep(0, num_trial)
  vec_pval4 = rep(0, num_trial)
  for(j in 1:num_trial) {
    # generates data
    mu1 = c(0, 0)
    mu2 = c(delta, 0)
    mu3 = c(2 * delta, 0)
    
    X1 = mvrnorm(n = n_each, mu1, diag(ss_diag))
    X2 = mvrnorm(n = n_each, mu2, diag(ss_diag))
    X3 = mvrnorm(n = n_each, mu3, diag(ss_diag))
    
    X = rbind(X1, X2, X3)
    
    # runs clustering
    hc = hclust(dist(X) ** 2, method = "average")
    cl = stats::cutree(hc, K) 
    
    # checks if the pair k1, k2 is under the alternative
    ind_pair = c(which(cl == k1), which(cl == k2))
    if (length(unique(cl_true[ind_pair])) == 1) {
      vec_pval1[j] = NA
      vec_pval2[j] = NA
      vec_pval3[j] = NA
      vec_pval4[j] = NA
    } else {
      # computes p-values
      vec_pval1[j] = fun_proposed_approx(X, K, k1, k2, ndraws, alpha)[1]
      ss_hat_all = fun_ss_hat_all(X, n, q)
      vec_pval2[j] = test_hier_clusters_exact(X, "average", hc, K, k1, k2, 
                                              sig = sqrt(ss_hat_all))$pval
      ss_hat_clustered = fun_ss_hat_clustered(X, n, q, cl, K)
      vec_pval3[j] = test_hier_clusters_exact(X, "average", hc, K, k1, k2, 
                                              sig = sqrt(ss_hat_clustered))$pval
      vec_pval4[j] = test_hier_clusters_exact(X, "average", hc, K, k1, k2, 
                                              sig = sqrt(ss))$pval
    }
  }
  # computes power
  tmp1 = fun_power(vec_pval1, thresh)
  tmp2 = fun_power(vec_pval2, thresh)
  tmp3 = fun_power(vec_pval3, thresh)
  tmp4 = fun_power(vec_pval4, thresh)
  mat_power[1, i] = tmp1[1]
  mat_power_eb[1, i] = tmp1[2]
  mat_power[2, i] = tmp2[1]
  mat_power_eb[2, i] = tmp2[2]
  mat_power[3, i] = tmp3[1]
  mat_power_eb[3, i] = tmp3[2]
  mat_power[4, i] = tmp4[1]
  mat_power_eb[4, i] = tmp4[2]
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
  theme(legend.position = "bottom",
        legend.direction = "vertical", 
        legend.text.align = 0) + 
  scale_color_manual(labels = c(TeX(r'(Proposed method with unknown $\sigma$)'),
                                TeX(r'(Gao et al.'s method with $\hat{\sigma}_{all}$)'), 
                                TeX(r'(Gao et al.'s method with $\hat{\sigma}_{clustered}$)'),
                                TeX(r'(Oracle $($Gao et al.'s method with true $\sigma)$)')), 
                     values = c("proposed" = "#C8375C",
                                "all" = "#3E9FB3", 
                                "clustered" = "#FFC918",
                                "oracle" = "#605D58")) + 
  ylim(0, 1) + 
  plot_details

##------------------------------------------------------------------------------
## Saves workspace

to_save = "workspaces/fig_power_K3_horizontal.RData"
if (dir.exists("workspaces") == 0) {
  dir.create("workspaces")
  save.image(file = to_save)
} else {
  save.image(file = to_save)
}

##------------------------------------------------------------------------------
## Saves plot

plot_comb = plot_comb = plot_sim_low | plot_sim_high | plot_power
to_save = paste0("plots/fig_power_K3_horizontal.png")
if (dir.exists("plots") == 0) {
  dir.create("plots")
  ggsave(plot_comb, file = to_save, width = 26, height = 9.2)
} else {
  ggsave(plot_comb, file = to_save, width = 26, height = 9.2)
}

##------------------------------------------------------------------------------