##------------------------------------------------------------------------------
## Imports libraries 

library(MASS)
library(intervals)
library(truncnorm)

##------------------------------------------------------------------------------
## Plot specifications

TITLE_SIZE = 30
AXIS_TITLE_SIZE = 30
AXIS_TEXT_SIZE = 30
LEGEND_TITLE_SIZE = 30
LEGENT_TEXT_SIZE = 30

plot_details = theme(plot.title = element_text(hjust = 0.5, 
                                               size = TITLE_SIZE),
                     axis.title = element_text(size = AXIS_TITLE_SIZE, 
                                               color = "black"),
                     axis.text = element_text(size = AXIS_TEXT_SIZE, 
                                              color = "black"),
                     axis.line = element_line(color = "black"), 
                     panel.background = element_rect(fill = 'white', 
                                                     color = 'white'),
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     legend.background = element_rect(fill = "transparent"), 
                     legend.title = element_text(size = LEGEND_TITLE_SIZE),
                     legend.text = element_text(size = LEGENT_TEXT_SIZE), 
                     legend.key = element_rect(color = NA, fill = NA, 
                                               size = unit(1, "cm")))

##------------------------------------------------------------------------------
## Functions for estimating variances

fun_ss_hat_clustered = function(X, n, q, cl, K) {
  # input(s): 
  #  - X: data matrix of dimensions n by q 
  #  - n: number of rows of X 
  #  - q: number of columns of X
  #  - cl: vector of cluster assignments, e.g. c(1, 1, 2, 3, 3, 2)
  #  - K: number of clusters produced by the clustering algorithm
  # output(s): 
  #  - variance estimator used for Gao et al. with sigmahat_clustered
  numer = 0 
  for (i in 1:K) {
    n_i = sum(cl == i)
    X_i = matrix(X[cl == i, ], ncol = q)
    center_i = colMeans(X_i)
    mat_center_i = matrix(rep(center_i, n_i), n_i, q, byrow = TRUE)
    numer = numer + sum((X[cl == i, ] - mat_center_i) ** 2)
  }
  to_return = numer / ((n - K) * q) 
  return(to_return)
}

fun_ss_hat_all = function(X, n, q) {
  # input(s): 
  #  - X: data matrix of dimensions n by q 
  #  - n: number of rows of X 
  #  - q: number of columns of X
  # output(s): 
  #  - Variance estimator used for Gao et al. with sigmahat_all
  center = colMeans(X)
  mat_center = matrix(rep(center, n), n, q, byrow = TRUE) 
  to_return = sum((X - mat_center) ** 2) / ((n - 1) * q)
  return(to_return)
}

##------------------------------------------------------------------------------
## Functions for decomposing X

fun_v = function(cl, k1, k2) {
  # input(s):
  # - cl: vector of cluster assignments, e.g. c(1, 1, 2, 3, 3, 2)
  # - k1: cluster index of one of the clusters we are testing
  # - k2: cluster index of the other cluster we are testing
  # output(s):
  # - v, which is a vector of length n
  n = length(cl)
  tmp1 = (cl == k1) / sum(cl == k1)
  tmp2 = (cl == k2) / sum(cl == k2)
  to_return = tmp1 - tmp2
  return(to_return)
}

fun_w = function(cl, k1, k2) {
  # input(s):
  # - cl: vector of cluster assignments, e.g. c(1, 1, 2, 3, 3, 2)
  # - k1: cluster index of one of the clusters we are testing
  # - k2: cluster index of the other cluster we are testing
  # output(s):
  # - w, which is a vector of length n
  n = length(cl)
  n1 = sum(cl == k1)
  n2 = sum(cl == k2)
  w = matrix(rep(1, n), n, 1)
  idx_other = setdiff(seq(1, n), c(which(cl == k1), which(cl == k2)))
  w[idx_other] = 0
  w = w / (n1 + n2)
  return(w)
}

fun_P0 = function(cl, k1, k2) {
  # input(s):
  # - cl: vector of cluster assignments, e.g. c(1, 1, 2, 3, 3, 2)
  # - k1: cluster index of one of the clusters we are testing
  # - k2: cluster index of the other cluster we are testing
  # output(s):
  # - P0, which is a projection matrix of size n by n 
  v = fun_v(cl, k1, k2)
  P0 = v %*% t(v) / sum(v ** 2)
  return(P0)
}

fun_P2 = function(cl, k1, k2) {
  # input(s):
  # - cl: vector of cluster assignments, e.g. c(1, 1, 2, 3, 3, 2)
  # - k1: cluster index of one of the clusters we are testing
  # - k2: cluster index of the other cluster we are testing
  # output(s):
  # - P2, which is a projection matrix of size n by n 
  n = length(cl)
  w = fun_w(cl, k1, k2)
  idx_other = setdiff(seq(1, n), c(which(cl == k1), which(cl == k2)))
  P2 = w %*% t(w) / sum(w ** 2)
  for (i in 1:length(idx_other)) {
    evec = rep(0, n)
    evec[idx_other[i]] = 1
    P2 = P2 + evec %*% t(evec) 
  }
  return(P2)
}

fun_P1 = function(cl, k1, k2) {
  # input(s):
  # - cl: vector of cluster assignments, e.g. c(1, 1, 2, 3, 3, 2)
  # - k1: cluster index of one of the clusters we are testing
  # - k2: cluster index of the other cluster we are testing
  # output(s):
  # - P1, which is a projection matrix of size n by n 
  n = length(cl)
  P1 = diag(n) - fun_P0(cl, k1, k2) - fun_P2(cl, k1, k2)
  return(P1)
}

##------------------------------------------------------------------------------
## Functions for computing p-value exactly 

fun_gen_X = function(n, q, ss, delta = 0) {
  # input(s):
  # - n: number of rows of X 
  # - q: number of columns of X
  # - ss: true variance 
  # - delta: signal strength that differentiates two true clusters,  
  #          0 if null is true 
  # output(s):
  # - data matrix of dimensions n by q 
  n_each = round(n / 2)
  mu1 = c(delta, rep(0, q - 1))
  X1 = mvrnorm(n = n_each, mu1, ss * diag(q))
  X2 = matrix(rnorm(n_each * q, 0, sd = sqrt(ss)), n_each, q)
  X = rbind(X1, X2)
  return(X)
}

fun_ts = function(X, cl, k1, k2) {
  # input(s):
  # - X: data matrix of dimensions n by q 
  # - cl: vector of cluster assignments, e.g. c(1, 1, 2, 3, 3, 2)
  # - k1: cluster index of one of the clusters we are testing
  # - k2: cluster index of the other cluster we are testing
  # output(s):
  # - test statistic that follows Beta distribution
  n = length(cl)
  P0X = fun_P0(cl, k1, k2) %*% X 
  P1X = fun_P1(cl, k1, k2) %*% X
  ts = (sqrt(sum(P0X ** 2)) / sqrt(sum(P0X ** 2) + sum(P1X ** 2))) ** 2
  return(ts)
}

fun_P1X_norm = function(X, cl, k1, k2) {
  # input(s):
  # - X: data matrix of dimensions n by q 
  # - cl: vector of cluster assignments, e.g. c(1, 1, 2, 3, 3, 2)
  # - k1: cluster index of one of the clusters we are testing
  # - k2: cluster index of the other cluster we are testing
  # output(s):
  # - Frobenius norm of P1X
  P1X = fun_P1(cl, k1, k2) %*% X 
  to_return = sqrt(sum(P1X ** 2))
  return(to_return)
}

fun_F_to_chi2 = function(x, n1, n2) {
  # input(s):
  # x: support of the F distribution
  # n1: the first parameter of the F distribution 
  # n2: the second parameter of the F distribution 
  # output(s):
  # - the corresponding support of the chi squared distribution
  # Remark(s): 
  # Uses SFA approximation as described in  
  # Li, Baibing, and Elaine B. Martin. "An approximation to the F
  # distribution using the chi-square distribution." Computational 
  # statistics & data analysis 40, no. 1 (2002): 21-26.
  numRow = nrow(x)
  numCol = ncol(x)
  to_return = matrix(0, numRow, numCol)
  for (i in 1:numRow) {
    for (j in 1:numCol) {
      if (is.infinite(x[i, j]) == 1) {
        to_return[i, j] = Inf
      } else {
        # n1 is k, n2 is n
        tmp1 = 2 * n2 + (n1 * x[i, j]) / 3 + (n1 - 2)
        tmp2 = 2 * n2 + (4 * n1 * x[i, j]) / 3
        lambda =  tmp1 / tmp2
        to_return[i, j] = lambda * n1 * x[i, j]
      }
    }
  }
  return(to_return)
}

fun_pval_sig_unknown = function(X, cl, S, k1, k2) {
  # input(s):
  # - X: data matrix of dimensions n by q 
  # - cl: vector of cluster assignments, e.g. c(1, 1, 2, 3, 3, 2)
  # - S: the truncation set
  # - k1: cluster index of one of the clusters we are testing
  # - k2: cluster index of the other cluster we are testing
  # output(s):
  # - p-value 
  n = length(cl)
  q = ncol(X)
  P1X_norm = fun_P1X_norm(X, cl, k1, k2)
  S_ = S / (P1X_norm * clusterpval:::norm_vec(fun_v(cl, k1, k2)))
  ts = fun_ts(X, cl, k1, k2) # test statistic that follows Beta distribution
  ts = (n - 2) * ts / (1 - ts) # converts to F statistic
  set_ts = Intervals(c(ts, Inf))
  set_denom = (n - 2) * S_ ** 2
  set_numer = interval_intersection(set_ts, set_denom)
  set_denom = fun_F_to_chi2(set_denom, q, (n - 2) * q)
  set_numer = fun_F_to_chi2(set_numer, q, (n - 2) * q)
  pval = clusterpval:::TChisqRatioApprox(q, set_numer, set_denom)
  return(pval)
}

fun_proposed_exact = function(X, k1 = 1, k2 = 2) {
  # input(s):
  # - X: data matrix of dimensions n by q 
  # - k1: cluster index of one of the clusters we are testing
  # - k2: cluster index of the other cluster we are testing
  # output(s):
  # - p-value
  hc = hclust(dist(X) ** 2, method = "average")
  cl = cutree(hc, 2) 
  S = clusterpval:::compute_S_average(X, hc, 2, k1, k2, dist(X) ** 2)
  pval = fun_pval_sig_unknown(X, cl, S, k1, k2)
  return(pval)
}

##------------------------------------------------------------------------------
## Functions for computing p-value approximately

fun_ratio_target_to_proposal = function(x, ts, m, q, alpha) {
  # input(s):
  # - x: value at which we're computing the density function
  # - ts: test statistic following the F distribution
  # - m: number of observations in the two clusters we are testing
  # - q: number of columns of X
  # - alpha: the sd of the proposal distribution
  # output(s):
  # - difference between the logs of the pdf value of the target distribution 
  #   and the pdf value of the proposal distribution 
  numer_f = dbeta(x, q / 2, (m * q - 2 * q) / 2, log = TRUE)
  denom_g = log(dtruncnorm(x, a = 0, b = 1, mean = ts, sd = alpha))
  to_return = numer_f - denom_g
  return(to_return)
}

fun_proposed_approx = function(X, K, k1, k2, ndraws, alpha) {
  # input(s):
  # - X: data matrix of dimensions n by q 
  # - k: number of clusters the clustering algorithm produces
  # - k1: cluster index of one of the clusters we are testing
  # - k2: cluster index of the other cluster we are testing
  # - ndraws: number of samples drawn in the importance sampling algorithm 
  # - alpha: the sd of the proposal distribution
  # output(s):
  # - p-value
  # - acceptance rate 
  hc = hclust(dist(X) ** 2, method = "average")
  cl = cutree(hc, K) 
  # computes useful quantities
  n = nrow(X)
  q = ncol(X)
  n1 = sum(cl == k1)
  n2 = sum(cl == k2)
  m = n1 + n2
  w = fun_w(cl, k1, k2)
  P0 = fun_P0(cl, k1, k2)
  P1 = fun_P1(cl, k1, k2)
  P2 = fun_P2(cl, k1, k2)
  # constructs test statistic
  psi0 = clusterpval:::norm_vec(P0 %*% X)
  psi1 = clusterpval:::norm_vec(P1 %*% X)
  ts = (psi0 / sqrt(psi0 ** 2 + psi1 ** 2)) ** 2 # follows Beta distribution
  # samples from the proposal distribution 
  samp = rtruncnorm(ndraws, a = 0, b = 1, mean = ts, sd = alpha) 
  # computes weights
  tmp_w = samp
  w = sapply(tmp_w, fun_ratio_target_to_proposal, ts = ts, m = m, 
             q = q, alpha = alpha)
  # computes part of h1 
  vec_h1 = rep(0, ndraws)
  vec_h1[samp >= ts] = 1
  # computes h2  
  t1 = sqrt(psi0 ** 2 + psi1 ** 2) 
  t2 = P0 %*% X / psi0 
  t3 = P1 %*% X / psi1
  vec_h2 = rep(0, ndraws)
  for (j in 1:ndraws) {
    samp_cur = samp[j]
    y = t1 * sqrt(samp_cur) * t2 + t1 * sqrt(1 - samp_cur) * t3 + P2 %*% X
    hc_y = hclust(dist(y) ** 2, method = "average")
    cl_y = cutree(hc_y, K) 
    if(clusterpval:::same_cl(cl, cl_y, K)) {
      vec_h2[j] = 1
    }
  }
  # computes the rest of h1
  vec_h1[vec_h2 == 0] = 0
  # computes the p-value
  w_new = w[vec_h2 == 1]
  w_new = w_new - max(w_new)
  pval = sum(exp(w_new[vec_h1[vec_h2 == 1] == 1]) / sum(exp(w_new)))
  acc_rate = sum(vec_h2) / ndraws
  return(c(pval, acc_rate))
}

##------------------------------------------------------------------------------
## Functions for computing empirical power

fun_power = function(vec_p, thresh) {
  # input(s):
  # - vec_p: vector of p-values of length num_trial 
  # - thresh: the level at which we are testing the null hypothesis 
  # output(s):
  # - empirical power
  # - standard error
  num_p = length(vec_p)
  indic = rep(0, num_p)
  indic[vec_p < thresh] = 1
  po = sum(indic) / num_p
  eb = sd(indic) / sqrt(num_p) 
  to_return = c(po, eb)
  return(to_return)
}

##------------------------------------------------------------------------------