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
ss_diag = ss * c(1, 1)
n = 120 
n_each = round(n / 3)
q = 2
delta = 10

# settings for clustering 
K = 3
k1 = 1
k2 = 2

##------------------------------------------------------------------------------
## Generates data

set.seed(2)

# generates original data
mu1 = c((delta / 2), (delta / 2))
mu2 = c(0, delta)
mu3 = c(delta, delta)
X1 = mvrnorm(n = n_each, mu1, diag(ss_diag))
X2 = mvrnorm(n = n_each, mu2, diag(ss_diag))
X3 = mvrnorm(n = n_each, mu3, diag(ss_diag))
X = rbind(X1, X2, X3)

# runs clustering on original data
hc = hclust(dist(X) ** 2, method = "average")
cl = cutree(hc, K) 

# saves original data and cluster assignments
df_orig = data.frame(x = X[, 1],
                     y = X[, 2],
                     col_vec = as.factor(cl))

# computes cluster centers
m = colMeans(X[which(cl == 1 | cl == 2), ])
m1 = colMeans(X[which(cl == 1), ])
m2 = colMeans(X[which(cl == 2), ])
m3 = colMeans(X[which(cl == 3), ])

# decomposes data into orthogonal components
P0 = fun_P0(cl, k1, k2)
P1 = fun_P1(cl, k1, k2)
P2 = fun_P2(cl, k1, k2)
P0X = P0 %*% X
P1X = P1 %*% X
P2X = P2 %*% X

# generates new data 1: 
# original data with first component scaled by 2
X1 = 2 * P0X + P1X + P2X

# computes cluster centers of new data 1
m_1 = colMeans(X1[which(cl == 1 | cl == 2), ])
m1_1 = colMeans(X1[which(cl == 1), ])
m2_1 = colMeans(X1[which(cl == 2), ])
m3_1 = colMeans(X1[which(cl == 3), ])

# saves new data 1 and cluster assignments
df1 = data.frame(x = X1[, 1],
                 y = X1[, 2],
                 col_vec = as.factor(cl))

# generates new data 2: 
# original data with second component scaled by 2
X2 = P0X + 2 * P1X + P2X

# computes cluster centers of new data 2
m_2 = colMeans(X2[which(cl == 1 | cl == 2), ])
m1_2 = colMeans(X2[which(cl == 1), ])
m2_2 = colMeans(X2[which(cl == 2), ])
m3_2 = colMeans(X2[which(cl == 3), ])

# saves new data 2 and cluster assignments
df2 = data.frame(x = X2[, 1],
                 y = X2[, 2],
                 col_vec = as.factor(cl))

##------------------------------------------------------------------------------
## Plots data 

# plots original data
plot_orig = ggplot(df_orig, aes(x = x, y = y)) + 
  geom_point(aes(color = col_vec), size = 3) + 
  geom_point(aes(x = m[1], y = m[2]), color = "blue", size = 3, 
             shape = 17) + 
  geom_point(aes(x = m1[1], y = m1[2]), color = "black", size = 3) + 
  geom_point(aes(x = m2[1], y = m2[2]), color = "black", size = 3) + 
  geom_point(aes(x = m3[1], y = m3[2]), color = "black", size = 3) + 
  scale_color_manual(values = c("#C8375C", "#3E9FB3", "#605D58")) + 
  labs(title = TeX(r"(\textit{$P_0X+P_1X+P_2X$})"),
       color = NULL, 
       x = NULL, y = NULL) +
  theme(legend.position = "none") + 
  xlim(c(-5, 15)) + 
  ylim(c(-1, 15)) + 
  plot_details

# plots new data 1
plot1 = ggplot(df1, aes(x = x, y = y)) + 
  geom_point(aes(color = col_vec), size = 3) + 
  geom_point(aes(x = m_1[1], y = m_1[2]), color = "blue", size = 3,
             shape = 17) + 
  geom_point(aes(x = m1_1[1], y = m1_1[2]), color = "black", size = 3) + 
  geom_point(aes(x = m2_1[1], y = m2_1[2]), color = "black", size = 3) + 
  geom_point(aes(x = m3_1[1], y = m3_1[2]), color = "black", size = 3) + 
  scale_color_manual(values = c("#C8375C", "#3E9FB3", "#605D58")) + 
  labs(title = TeX(r"(\textit{$2P_0X+P_1X+P_2X$})"), 
       color = NULL, 
       x = NULL, y = NULL) +
  theme(legend.position = "none") + 
  xlim(c(-5, 15)) + 
  ylim(c(-1, 15)) + 
  plot_details

# plots new data 2
plot2 = ggplot(df2, aes(x = x, y = y)) + 
  geom_point(aes(color = col_vec), size = 3) + 
  geom_point(aes(x = m_2[1], y = m_2[2]), color = "blue", size = 3,
             shape = 17) + 
  geom_point(aes(x = m1_2[1], y = m1_2[2]), color = "black", size = 3) + 
  geom_point(aes(x = m2_2[1], y = m2_2[2]), color = "black", size = 3) + 
  geom_point(aes(x = m3_2[1], y = m3_2[2]), color = "black", size = 3) + 
  scale_color_manual(values = c("#C8375C", "#3E9FB3", "#605D58")) + 
  labs(title =  TeX(r"(\textit{$P_0X+2P_1X+P_2X$})"), 
       color = NULL, 
       x = NULL, y = NULL) +
  theme(legend.position = "none") + 
  xlim(c(-5, 15)) + 
  ylim(c(-1, 15)) + 
  plot_details

##------------------------------------------------------------------------------
## Saves workspace

to_save = "workspaces/fig_decomposition.RData"
if (dir.exists("workspaces") == 0) {
  dir.create("workspaces")
  save.image(file = to_save)
} else {
  save.image(file = to_save)
}

##------------------------------------------------------------------------------
## Saves plot

plot_comb = (plot_orig | plot1 | plot2)
to_save = "plots/fig_decomposition.png"
if (dir.exists("plots") == 0) {
  dir.create("plots")
  ggsave(plot_comb, file = to_save, width = 25.5, height = 6.5)
} else {
  ggsave(plot_comb, file = to_save, width = 25.5, height = 6.5)
}

##------------------------------------------------------------------------------