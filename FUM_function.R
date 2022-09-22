FCM <- function(x, K, mybeta = 2, nstart = 1, iter_max = 100, eps = 1e-06) {

  ## modified time: 2018-03-23

  FCM_onetime <- function(x, init_centers, mybeta = 2, iter_max = 100, eps = 1e-06) {
    n = dim(x)[1]
    d = dim(x)[2]
    g = init_centers
    K = dim(g)[1]
    histJ = c()
    pasfini = 1
    Jold = Inf
    D = matrix(0, n, K)
    for (j in 1:K) {
      D[, j] = rowSums(sweep(x, 2, g[j, ], "-")^2)
    }
    iter = 1
    J_old = Inf
    while (pasfini) {
      s = (1/(D + eps))^(1/(mybeta - 1))
      u = s/(s %*% matrix(1, K, K))
      t1 = t(u^mybeta) %*% x
      t2 = t(u^mybeta) %*% matrix(1, n, d)
      V = t1/t2
      g = V
      D = matrix(0, n, K)
      for (j in 1:K) {
        D[, j] = rowSums(sweep(x, 2, g[j, ], "-")^2)
      }
      J = sum(u^mybeta * D)
      pasfini = abs(J - Jold) > 0.001 && (iter < iter_max)
      Jold = J
      histJ = c(histJ, J)
      iter = iter + 1
    }
    cluster_id = apply(u, 1, which.max)
    re = list(u, J, histJ, g, cluster_id)
    names(re) = c("u", "J", "histJ", "g", "cluster_id")
    return(re)
  }
  x = as.matrix(x)
  seeds = 1:nrow(x)
  id = sample(seeds, K)
  g = as.matrix(x[id, ])
  re_best = FCM_onetime(x = x, init_centers = g, mybeta = mybeta, iter_max = iter_max, eps = eps)
  if (nstart > 1) {
    minJ = 0
    i = 2
    while (i <= nstart) {
      init_centers_id = sample(seeds, K)
      init_centers = as.matrix(x[init_centers_id, ])
      run = FCM_onetime(x, init_centers = init_centers, mybeta = mybeta, iter_max = iter_max)
      if (run$J <= re_best$J) {
        re_best = run
      }
      i = i + 1
    }
  }
  return(re_best)
} 