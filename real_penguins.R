##------------------------------------------------------------------------------
## Imports libraries 

library(fastcluster)
library(clusterpval)
library(ggplot2)
library(patchwork)
library(xtable)
library(palmerpenguins)
source("functions.R")

##------------------------------------------------------------------------------
## Specifies settings 

# settings for clustering 
K = 6

# settings for importance sampling
ndraws = 8000 # number of samples drawn
alpha = 0.05 # sd of the proposal distribution

##------------------------------------------------------------------------------
## Subsets the penguins dataset to female penguins and to 
## flipper and bill lengths

penguins = na.omit(penguins)
penguins = penguins[penguins$sex == "female", ]
X_pre = as.matrix(penguins[, c(3, 5)])
true_y = as.vector(as.matrix(penguins[, 1]))
n = dim(X_pre)[1]
q = dim(X_pre)[2]

##------------------------------------------------------------------------------
## Visualizes original data with cluster assignments generated with 
## standardized data

X = apply(X_pre, 2, function (x) {(x - mean(x)) / sd(x)})

# runs clustering
hc = hclust(dist(X) ** 2, method = "average")
cl = stats::cutree(hc, K) 

# visualizes data 
df = data.frame(x = X_pre[, 2],
                y = X_pre[, 1],
                col_vec = as.factor(cl),
                shape_vec = as.factor(rep(c("Adelie", "Gentoo", "Chinstrap"), 
                                          each = round(n / 3))))

plot_penguins = ggplot(df, aes(x = x, y = y)) + 
  geom_point(aes(shape = shape_vec, color = col_vec), size = 5) + 
  scale_shape_manual(values = c("Adelie" = 15, "Gentoo" = 16, 
                                "Chinstrap" = 17)) + 
  scale_color_manual(values = c("1" = "#898989", 
                                "2" = "#363C44", 
                                "3" = "#006992", 
                                "4" = "#A1A0D2", 
                                "5" = "#64379E", 
                                "6" = "#67B3C9")) + 
  labs(title = "Visualization of Data", 
       color = "Clusters Formed", shape = "True Clusters", 
       x = "Flipper length (mm)", y = "Bill length (mm)") +
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 20),
        axis.title = element_text(size = 20, 
                                  color = "black"),
        axis.text = element_text(size = 20, 
                                 color = "black"),
        axis.line = element_line(color = "black"), 
        panel.background = element_rect(fill = 'white', 
                                        color = 'white'),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = "transparent"), 
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20), 
        legend.key = element_rect(color = NA, fill = NA, 
                                  size = unit(1, "cm")))
plot_penguins

##------------------------------------------------------------------------------
# Computes p-values 

# sets cluster pairs to test
vec_ind1 = c(1, 1, 4)
vec_ind2 = c(2, 5, 5)
l = length(vec_ind1)

mat_sim = matrix(0, 3, l)

# computes p-values using proposed method with unknown sigma 
set.seed(2)
count = 1
for (i in 1:l) {
  k1 = vec_ind1[i]
  k2 = vec_ind2[i]
  mat_sim[1, i] = fun_proposed_approx(X, K, k1, k2, ndraws, alpha)[1]
}

# computes p-values using Gao et al. with sigmahat_all
set.seed(2)
# estimates variance
ss_hat_all = fun_ss_hat_all(X, n, q)
count = 1
for (i in 1:l) {
  k1 = vec_ind1[i]
  k2 = vec_ind2[i]
  # computes the p-value 
  mat_sim[2, i] = test_hier_clusters_exact(X, "average", hc, K, k1, k2, 
                                           sig = sqrt(ss_hat_all))$pval
}

# computes p-values using Gao et al. with sigmahat_clustered
set.seed(2)
# estimates variance
ss_hat_clustered = fun_ss_hat_clustered(X, n, q, cl, K)
count = 1
for (i in 1:l) {
  k1 = vec_ind1[i]
  k2 = vec_ind2[i]
  # computes the p-value 
  mat_sim[3, i] = test_hier_clusters_exact(X, "average", hc, K, k1, k2, 
                                           sig = sqrt(ss_hat_clustered))$pval
}

##------------------------------------------------------------------------------
## Prints results
## Citations: the codes for creating a LATEX table are from the 
##            xtable documentation.

# prints results
xt_sim = xtable(data.frame(c1 = mat_sim[, 1], 
                           c2 = mat_sim[, 2], 
                           c3 = mat_sim[, 3],
                           row.names = c("proposed", 
                                         "all", 
                                         "clustered")), 
                display = c("s", "g", "g", "g"))
names(xt_sim) = c("(1,2)", "(1,5)", "(4,5)")
print(xt_sim, math.style.exponents = TRUE)

##------------------------------------------------------------------------------
## Saves workspace

to_save = "workspaces/real_penguins.RData"
if (dir.exists("workspaces") == 0) {
  dir.create("workspaces")
  save.image(file = to_save)
} else {
  save.image(file = to_save)
}

##------------------------------------------------------------------------------
## Saves plot

to_save = "plots/real_penguins.png"
if (dir.exists("plots") == 0) {
  dir.create("plots")
  ggsave(plot_penguins, file = to_save, width = 10, height = 6)
} else {
  ggsave(plot_penguins, file = to_save, width = 10, height = 6)
}

##------------------------------------------------------------------------------
## Saves xtable result

to_save = "results/real_penguins.txt"
if (dir.exists("results") == 0) {
  dir.create("results")
  print(xt_sim, file = to_save)
} else {
  print(xt_sim, file = to_save)
}

##------------------------------------------------------------------------------