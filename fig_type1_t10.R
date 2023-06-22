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
n = 30
vec_q = c(2, 2)
df_t = 10

# settings for clustering 
vec_K = c(2, 3)
k1 = 1
k2 = 2

# settings for importance sampling
ndraws = 8000 # number of samples drawn
alpha = 0.05 # sd of the proposal distribution

# other settings
num_trial = 2000 # number of trials

##------------------------------------------------------------------------------
## Creates a data set consistent with the null for first setting

K = 2
q = 2

set.seed(2)
X = matrix(rt(n * q, df_t), n, q)

# runs clustering
hc = hclust(dist(X) ** 2, method = "average")
cl = stats::cutree(hc, K)

# visualizes
df_sim_1 = data.frame(x = X[, 1],
                      y = X[, 2],
                      col_vec = as.factor(cl))
plot_sim_1 = ggplot(df_sim_1, aes(x = x, y = y)) +
  geom_point(aes(color = col_vec), size = 5, stroke = 3) +
  scale_color_manual(values = c("1" = "#898989",
                                "2" = "#64379E")) +
  labs(title = TeX(r"(Visualization for \textit{$K=2$})"),
       color = element_blank(), shape = "True Clusters",
       x = NULL, y = NULL) +
  theme(legend.position = "none") + 
  plot_details

##------------------------------------------------------------------------------
## Creates a data set consistent with the null for second setting

K = 3
q = 2

set.seed(2)
X = matrix(rt(n * q, df_t), n, q)

# runs clustering
hc = hclust(dist(X) ** 2, method = "average")
cl = stats::cutree(hc, K)

# visualizes
df_sim_2 = data.frame(x = X[, 1],
                      y = X[, 2],
                      col_vec = as.factor(cl))
plot_sim_2 = ggplot(df_sim_2, aes(x = x, y = y)) +
  geom_point(aes(color = col_vec), size = 5, stroke = 3) +
  scale_color_manual(values = c("1" = "#898989",
                                "2" = "#64379E",
                                "3" = "#67B3C9")) +
  labs(title = TeX(r"(Visualization for \textit{$K=3$})"),
       color = element_blank(), shape = "True Clusters",
       x = NULL, y = NULL) +
  theme(legend.position = "none") + 
  plot_details

##------------------------------------------------------------------------------
## Creates QQ plots of p-values

list_mat = list()
for (l in 1:2) {
  K = vec_K[l]
  q = vec_q[l]
  mat_pval = matrix(0, 3, num_trial)
  if (K == 2) {
    set.seed(2)
    for (j in 1:num_trial) {
      X = matrix(rt(n * q, df_t), n, q)
      mat_pval[1, j] = fun_proposed_exact(X)
    }
  } else {
    set.seed(2)
    alpha = 0.02
    for (j in 1:num_trial) {
      X = matrix(rt(n * q, df_t), n, q)
      mat_pval[1, j] = fun_proposed_approx(X, K, k1, k2, ndraws, alpha)[1]
    }
  }
  set.seed(2)
  for (j in 1:num_trial) {
    X = matrix(rt(n * q, df_t), n, q)
    # runs clustering
    hc = hclust(dist(X) ** 2, method = "average")
    cl = cutree(hc, K)
    # estimates variance 
    ss_hat_all = fun_ss_hat_all(X, n, q)
    ss_hat_clustered = fun_ss_hat_clustered(X, n, q, cl, K)
    mat_pval[2, j] = test_hier_clusters_exact(X, "average", hc, K, k1, k2, 
                                              sig = sqrt(ss_hat_all))$pval
    mat_pval[3, j] = test_hier_clusters_exact(X, "average", hc, K, k1, k2, 
                                              sig = sqrt(ss_hat_clustered))$pval
  }
  list_mat[[l]] = mat_pval
}

##------------------------------------------------------------------------------
## Creates plots

df1 = data.frame(pval = c(t(list_mat[[1]])), 
                 type = as.factor(rep(c("proposed", "all", 
                                        "clustered"), 
                                      each = num_trial)))

plot1 = ggplot(df1, aes(sample = pval, color = type)) + 
  stat_qq(distribution = qunif) + 
  geom_abline(intercept = 0, slope = 1, col = "black", 
              size = 1, linetype = "dashed") + 
  labs(title = TeX(r"(QQ Plot for \textit{$K=2$})"),
       x = "Theoretical Quantiles", y = "Empirical Quantiles", 
       color = element_blank()) + 
  xlim(c(0, 1)) + 
  ylim(c(0, 1)) + 
  scale_color_manual(labels = c(TeX(r'(Proposed method with unknown $\sigma$)'),
                                TeX(r'(Gao et al.'s method with $\hat{\sigma}_{all}$)'), 
                                TeX(r'(Gao et al.'s method with $\hat{\sigma}_{clustered}$)')),
                     values = c("proposed" = "#C8375C",
                                "all" = "#3E9FB3", 
                                "clustered" = "#FFC918")) + 
  theme(legend.position = "none", 
        legend.text.align = 0) + 
  guides(color = guide_legend(override.aes = list(size = 3))) + 
  plot_details

df2 = data.frame(pval = c(t(list_mat[[2]])), 
                 type = as.factor(rep(c("proposed", "all", 
                                        "clustered"), 
                                      each = num_trial)))

plot2 = ggplot(df2, aes(sample = pval, color = type)) + 
  stat_qq(distribution = qunif) + 
  geom_abline(intercept = 0, slope = 1, col = "black", 
              size = 1, linetype = "dashed") + 
  labs(title = TeX(r"(QQ Plot for \textit{$K=3$})"),
       x = "Theoretical Quantiles", y = "Empirical Quantiles", 
       color = element_blank()) + 
  xlim(c(0, 1)) + 
  ylim(c(0, 1)) + 
  scale_color_manual(labels = c(TeX(r'(Proposed method with unknown $\sigma$)'),
                                TeX(r'(Gao et al.'s method with $\hat{\sigma}_{all}$)'), 
                                TeX(r'(Gao et al.'s method with $\hat{\sigma}_{clustered}$)')),
                     values = c("proposed" = "#C8375C",
                                "all" = "#3E9FB3", 
                                "clustered" = "#FFC918")) + 
  theme(legend.position = "none", 
        legend.text.align = 0) + 
  guides(color = guide_legend(override.aes = list(size = 3))) + 
  plot_details

##------------------------------------------------------------------------------
## Saves workspace

to_save = "workspaces/fig_type1_t10.RData"
if (dir.exists("workspaces") == 0) {
  dir.create("workspaces")
  save.image(file = to_save)
} else {
  save.image(file = to_save)
}

##------------------------------------------------------------------------------
## Saves plot

plot_comb1 = plot_sim_1 / plot_sim_2
plot_comb2 = (plot1 / plot2) & 
  theme(legend.position = "right",
        legend.justification = "right",
        legend.direction = "vertical") &
  guides(color = guide_legend(override.aes = list(size = 3))) 
plot_comb = plot_comb1 | plot_comb2
to_save = "plots/fig_type1_t10.png"
if (dir.exists("plots") == 0) {
  dir.create("plots")
  ggsave(plot_comb + plot_layout(guides = "collect"), 
         file = to_save, width = 27, height = 13)
} else {
  ggsave(plot_comb + plot_layout(guides = "collect"), 
         file = to_save, width = 27, height = 13)
}

##------------------------------------------------------------------------------