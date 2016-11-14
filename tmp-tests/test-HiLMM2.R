require(HiLMM)

data_sim <- data_simu(n = 1e3, N = 1e4, eta_star = 0.98, q = 0.05)
W <- data_sim$W
Y <- data_sim$Y

n <- nrow(W)
N <- ncol(W)

Z <- scale(W, center = TRUE, scale = TRUE)
M <- tcrossprod(Z) / N
eigs <- eigen(M, symmetric = TRUE)
lambda <- eigs$values[-n]
Y_tilde <- crossprod(eigs$vectors[, -n], Y)

tmp1 <- lambda - 1
Y_tilde2 <- Y_tilde^2

tosolve <- function(eta) {
  tmp2 <- (eta * tmp1 + 1)
  tmp3 <- tmp1 / tmp2
  tmp4 <- Y_tilde2 / tmp2

  tmp <- sum(tmp4 * tmp3) * n - sum(tmp4) * sum(tmp3)
  tmp^2
}

eta <- seq(0.001, 0.999, 0.001)
tmp2 <- eta %o% tmp1 + 1

test <- sapply(seq(0.01, 0.99, 0.01), tosolve)
plot(seq(0.01, 0.99, 0.01), test, type = "l")
test2 <- optimize(tosolve, lower = 0.01, upper = 0.99)
print(test2$minimum)
test3 <- estim_herit2(Y, W)
print(test3$heritability)
