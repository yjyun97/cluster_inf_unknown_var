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

# settings for clustering 
K = 2
k1 = 1
k2 = 2

# other settings
thresh = 0.05 # significance level
num_trial = 2000 # number of trials

##------------------------------------------------------------------------------
## Creates a data set consistent with the null

set.seed(6)
# generates data
X = fun_gen_X(n, q, ss)
# runs clustering
hc = hclust(dist(X) ** 2, method = "average")
cl = cutree(hc, K) 
# estimates variance
ss_hat_all = fun_ss_hat_all(X, n, q)
ss_hat_clustered = fun_ss_hat_clustered(X, n, q, cl, K)
# visualizes clusters
df = data.frame(x = X[, 1],
                y = X[, 2],
                col_vec = as.factor(cl),
                shape_vec = as.factor(rep(1, n)))
plot_null = ggplot(df, aes(x = x, y = y)) + 
  geom_point(aes(shape = shape_vec, color = col_vec), size = 5) + 
  scale_shape_manual(values = 15) + 
  scale_color_manual(values = c("1" = "#64379E",
                                "2" = "#898989")) + 
  labs(title = TeX(paste0("$\\hat{\\sigma}_{all}=$",
                          round(sqrt(ss_hat_all), 2),
                          " and ",
                          "$\\hat{\\sigma}_{clustered}=$",
                          round(sqrt(ss_hat_clustered), 2))),
       color = "Clusters Formed", shape = "True Clusters", 
       x = NULL, y = "Null") +
  theme(legend.position = "none") + 
  plot_details + 
  theme(axis.title = element_text(size = 40))
  
##------------------------------------------------------------------------------
## Creates a data set consistent with the alternative

set.seed(6)
# generates data
X = fun_gen_X(n, q, ss, 7)
# runs clustering
hc = hclust(dist(X) ** 2, method = "average")
cl = cutree(hc, K) 
# estimates variance 
ss_hat_all = fun_ss_hat_all(X, n, q)
ss_hat_clustered = fun_ss_hat_clustered(X, n, q, cl, K)
# visualizes clusters
df = data.frame(x = X[, 1],
                y = X[, 2],
                col_vec = as.factor(cl),
                shape_vec = as.factor(rep(c(1, 2), each = round(n / 2))))
plot_alter = ggplot(df, aes(x = x, y = y)) + 
  geom_point(aes(shape = shape_vec, color = col_vec), size = 5) + 
  scale_color_manual(values = c("1" = "#64379E",
                                "2" = "#898989")) + 
  scale_shape_manual(values = c(15, 16)) + 
  labs(title = TeX(paste0("$\\hat{\\sigma}_{all}=$",
                          round(sqrt(ss_hat_all), 2),
                          " and ",
                          "$\\hat{\\sigma}_{clustered}=$",
                          round(sqrt(ss_hat_clustered), 2))), 
       color = "Clusters Formed", shape = "True Clusters", 
       x = NULL, y = "Alternative") +
  theme(legend.position = "none") + 
  plot_details + 
  theme(axis.title = element_text(size = 40))

##------------------------------------------------------------------------------
## Computes p-values

vec_pval_all = rep(0, num_trial)
vec_pval_clustered = rep(0, num_trial)

set.seed(1)
for (i in 1:num_trial) {
  #generates data
  X = fun_gen_X(n, q, ss)
  # runs clustering
  hc = hclust(dist(X) ** 2, method = "average")
  cl = cutree(hc, K) 
  # estimates variance
  ss_hat_all = fun_ss_hat_all(X, n, q)
  ss_hat_clustered = fun_ss_hat_clustered(X, n, q, cl, K)
  # computes the p-value 
  vec_pval_all[i] = test_hier_clusters_exact(X, "average", hc, K, 1, 2, 
                                                sig = sqrt(ss_hat_all))$pval
  vec_pval_clustered[i] = test_hier_clusters_exact(X, "average", hc, K, 1, 2, 
                                                sig = sqrt(ss_hat_clustered))$pval
}

##------------------------------------------------------------------------------
## Creates QQ plot of p-values

df_pval = data.frame(all = vec_pval_all, 
                     clustered = vec_pval_clustered)

plot_pval_all = ggplot(df_pval, aes(sample = all)) + 
  stat_qq(distribution = qunif, color = "#3E9FB3",) + 
  geom_abline(slope = 1, intercept = 0, col = "black", 
              size = 1, linetype = "dashed") +
  labs(title = TeX(r'(QQ Plot of p-value with $\hat{\sigma}_{all}$)'), 
       x = "Theoretical Quantiles", y = "Empirical Quantiles") + 
  theme(legend.position = "none") + 
  xlim(c(0, 1)) + 
  ylim(c(0, 1)) + 
  plot_details

plot_pval_clustered = ggplot(df_pval, aes(sample = clustered)) + 
  stat_qq(distribution = qunif, color = "#FFC918") + 
  geom_abline(slope = 1, intercept = 0, col = "black", 
              size = 1, linetype = "dashed") +
  labs(title = TeX(r'(QQ Plot of p-value with $\hat{\sigma}_{clustered}$)'), 
       x = "Theoretical Quantiles", y = "Empirical Quantiles") + 
  theme(legend.position = "none") + 
  xlim(c(0, 1)) + 
  ylim(c(0, 1)) + 
  plot_details

##------------------------------------------------------------------------------
## Computes empirical power

mat_power = matrix(0, 3, ld) 
mat_power_eb = matrix(0, 3, ld)

for (i in 1:ld) {
  delta = vec_delta[i]
  set.seed(1)
  vec_pval1 = rep(0, num_trial)
  vec_pval2 = rep(0, num_trial)
  vec_pval3 = rep(0, num_trial)
  for(j in 1:num_trial) {
    # generates data
    X = fun_gen_X(n, q, ss, delta)
    # runs clustering
    hc = hclust(dist(X) ** 2, method = "average")
    cl = stats::cutree(hc, K) 
    ss_hat_all = fun_ss_hat_all(X, n, q)
    vec_pval1[j] = test_hier_clusters_exact(X, "average", hc, K, 1, 2,
                                            sig = sqrt(ss_hat_all))$pval
    ss_hat_clustered = fun_ss_hat_clustered(X, n, q, cl, K)
    vec_pval2[j] = test_hier_clusters_exact(X, "average", hc, K, 1, 2,
                                            sig = sqrt(ss_hat_clustered))$pval
    vec_pval3[j] = test_hier_clusters_exact(X, "average", hc, K, 1, 2,
                                            sig = sqrt(ss))$pval
  }
  tmp1 = fun_power(vec_pval1, thresh)
  tmp2 = fun_power(vec_pval2, thresh)
  tmp3 = fun_power(vec_pval3, thresh)
  mat_power[1, i] = tmp1[1]
  mat_power_eb[1, i] = tmp1[2]
  mat_power[2, i] = tmp2[1]
  mat_power_eb[2, i] = tmp2[2]
  mat_power[3, i] = tmp3[1]
  mat_power_eb[3, i] = tmp3[2]
}

##------------------------------------------------------------------------------
## Creates plot of empirical power 

stack = c(mat_power[1, ], mat_power[3, ])
stack_eb = c(mat_power_eb[1, ], mat_power_eb[3, ])

df_all = data.frame(delta = rep(vec_delta, 2),  
                power = stack, 
                case = factor(rep(c("all", "Oracle"), each = ld)), 
                leb = stack - stack_eb, 
                ueb = stack + stack_eb)

plot_power_all = ggplot(df_all, aes(delta, power, group = case)) + 
  geom_errorbar(aes(ymin = leb, ymax = ueb, col = case), 
                width = 0.3, size = 2) + 
  geom_line(aes(col = case), linetype = "solid", size = 2) + 
  geom_point(aes(col = case), size = 2) + 
  labs(title = TeX(r'(Empirical Power with $\hat{\sigma}_{all}$)'),
       x = TeX(r'($\delta$)'), y = "Empirical Power", color = NULL, 
       linetype = NULL) + 
  theme(legend.position = "bottom",
        legend.direction = "vertical", 
        legend.justification = "left",
        legend.text.align = 0) + 
  scale_color_manual(labels = c(TeX(r'(Gao et al.'s method with $\hat{\sigma}_{all}$)'), 
                                TeX(r'(Oracle $($Gao et al.'s method with true $\sigma)$)')),
                     values = c("all" = "#3E9FB3", 
                                "Oracle" = "#605D58")) + 
  ylim(0, 1) + 
  plot_details

stack = c(mat_power[2, ], mat_power[3, ])
stack_eb = c(mat_power_eb[2, ], mat_power_eb[3, ])

df_clustered = data.frame(delta = rep(vec_delta, 2),  
                       power = stack, 
                       case = factor(rep(c("clustered", "Oracle"), each = ld)), 
                       leb = stack - stack_eb, 
                       ueb = stack + stack_eb)

plot_power_clustered = ggplot(df_clustered, aes(delta, power, group = case)) + 
  geom_errorbar(aes(ymin = leb, ymax = ueb, col = case), 
                width = 0.3, size = 2) + 
  geom_line(aes(col = case), linetype = "solid", size = 2) + 
  geom_point(aes(col = case), size = 2) + 
  labs(title = TeX(r'(Empirical Power with $\hat{\sigma}_{clustered}$)'), 
       x = TeX(r'($\delta$)'), y = "Empirical Power", color = NULL, 
       linetype = NULL) + 
  theme(legend.position = "bottom",
        legend.direction = "vertical", 
        legend.justification = "left",
        legend.text.align = 0) + 
  scale_color_manual(labels = c(TeX(r'(Gao et al.'s method with $\hat{\sigma}_{clustered}$)'), 
                                TeX(r'(Oracle $($Gao et al.'s method with true $\sigma)$)')),
                       values = c("clustered" = "#FFC918", 
                                "Oracle" = "#605D58")) + 
  ylim(0, 1) + 
  plot_details

##------------------------------------------------------------------------------
## Saves workspace

to_save = "workspaces/sim_motivation.RData"
if (dir.exists("workspaces") == 0) {
  dir.create("workspaces")
  save.image(file = to_save)
} else {
  save.image(file = to_save)
}

##------------------------------------------------------------------------------
## Saves plot

plot_comb = (plot_null | plot_pval_clustered | plot_pval_all) / 
  (plot_alter | plot_power_clustered | plot_power_all)
to_save = "plots/sim_motivation.png"
if (dir.exists("plots") == 0) {
  dir.create("plots")
  ggsave(plot_comb, file = to_save, width = 28, height = 14)
} else {
  ggsave(plot_comb, file = to_save, width = 28, height = 14)
}

##------------------------------------------------------------------------------