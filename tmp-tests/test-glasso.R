# Adapted from https://github.com/JClavel/glassoFast/blob/master/src/glassofast.f90
Rcpp::sourceCpp("tmp-tests/test-glasso.cpp")

abs2 <- abs3 <- abs4 <- abs

glassofast2 <- function(S, L, maxit = 100, thr = 1e-4) {

  shr <- sum(abs2(S)) - sum(diag(S))
  n <- ncol(S)
  shr = thr*shr/(n-1)
  thrLasso = max(shr/n, 2e-16)

  WXj <- rep(0, n)
  W <- S + 0
  X <- S * 0

  Wd <- diag(S) + L
  diag(W) <- Wd

  for (iter in seq_len(maxit)) {

    print(iter)

    # dw <- 0
    # for (j in 1:n) {
    #   # WXj[] <- 0
    #   # for (i in 1:n) {
    #   #   if (X[i, j] != 0) {
    #   #     WXj = WXj + W[, i] * X[i, j]
    #   #   }
    #   # }
    #   ind <- which(X[, j] != 0)
    #   WXj <- W[, ind, drop = FALSE] %*% X[ind, j]
    #
    #   # repeat {
    #   #
    #   #   dlx = 0.0
    #   #   for (i in 1:n) {
    #   #     if (j != i) {
    #   #       a = S[i, j] - WXj[i] + Wd[i] * X[i, j]
    #   #       b = abs3(a) - L
    #   #       c = `if`(b > 0, sign(a) * b / Wd[i], 0)
    #   #       delta = c - X[i,j]
    #   #       if (delta != 0) {
    #   #         X[i,j] = c
    #   #         # WXj = WXj + W[, i] * delta
    #   #         update_WXj(WXj, W, i - 1L, delta)
    #   #         dlx = max(dlx, abs4(delta))
    #   #       }
    #   #     }
    #   #   }
    #   #
    #   #   if (dlx < thrLasso) break
    #   # }
    #   inner_loop(S, W, X, WXj, Wd, n, j - 1L, L, thrLasso)
    #
    #   WXj[j] = Wd[j]
    #   # dw = max(dw, sum(abs(WXj - W[, j])))
    #   dw = max(dw, abs_dist(WXj, W[, j]))
    #   W[, j] = W[j, ] <- WXj
    # }
    dw <- inner_loop(S, W, X, WXj, Wd, n, L, thrLasso)

    if (dw < shr) break
  }

  for (i in 1:n) {
    tmp <- 1 / drop(Wd[i] - crossprod(X[, i], W[, i]))
    X[, i] = -tmp * X[, i]
    X[i, i] = tmp
  }
  for (i in seq_len(n - 1)) {
    X[(i+1):n,i] = (X[(i+1):n,i] + X[i,(i+1):n])/2
    X[i,(i+1):n] = X[(i+1):n,i]
  }

  list(w = W, wi = X, niter = iter)
}

library(bigstatsr)
S <- cov(big_attachExtdata()[, 1:1000])

system.time(
  test <- glassofast2(S, L = 0.005)
) # 40 sec

hist(test$w - S)
mean(test$wi == 0)
# 0.050 -> 98.1%
# 0.010 -> 61.1% + weird distrib
# 0.005 -> 44.3% + weird distrib

system.time(
  test2 <- glassoFast::glassoFast(S, rho = 0.05, trace = TRUE)
) # 55 sec
str(test2)
all.equal(test, test2[-3])


microbenchmark::microbenchmark(
  {
    WXj[] <- 0
    for (i in 1:n) {
      if (X[i, j] != 0) {
        WXj = WXj + W[, i] * X[i, j]
      }
    }
    WXj
  },
  {
    ind <- which(X[, j] != 0)
    WXj <- W[, ind, drop = FALSE] %*% X[ind, j]
  },
  check = "equal"
)

S2 <- ifelse(S^2 > 5e-5, S, 0)
mean(S2 == 0)
system.time(
  test <- glassofast2(S2, L = 0.01)
) # 40 sec
hist(S2 - test$w)
norm(S2 - test$w, "F")

cov2cor(S2)[1:5, 1:5]
cov2cor(test$w)[1:5, 1:5]

RSpectra::eigs_sym(S2, k = 2, sigma = -1)$values
# -0.09737118 -0.09930095
RSpectra::eigs_sym(test$w, k = 2, sigma = -1)$values
# 0.044255635 0.009015971


system.time(
  test <- glassofast2(S, L = 0.001)
) # 405 sec
hist(S2 - test$w)
norm(S2 - test$w, "F")
mean(test$wi == 0) # 23%

system.time(
  test3 <- glassofast2(test$wi, L = 0.01)
) # 282 sec (stopped at 100 iter)

mean(test3$wi == 0) # 3.6%

hist(S - test3$wi)
norm(S - test3$wi, "F")

cov2cor(S)[1:5, 1:5]
cov2cor(test3$wi)[1:5, 1:5]

covg
