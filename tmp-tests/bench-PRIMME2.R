set.seed(1)
n <- 2e3
p <- 20e3
U <- matrix(0, n, 10); U[] <- rnorm(length(U))
V <- matrix(0, p, 10); V[] <- rnorm(length(V))
X <- tcrossprod(U, V) + 5 * matrix(rnorm(n * p), n, p)
system.time(svd0 <- svd(X, nu = 10))
system.time(svd1 <- RSpectra::svds(X, 10, opts = list(tol = 1e-4)))
system.time(svd2 <- PRIMME::svds(X, 10, isreal = TRUE))
plot(svd1$u, svd2$u, pch = 20)

log10(mean(abs(abs(svd0$u) - abs(svd1$u))))
log10(mean(abs(abs(svd0$u) - abs(svd2$u))))
