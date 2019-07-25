set.seed(1)
n <- 2e3
p <- 20e3
U <- matrix(0, n, 10); U[] <- rnorm(length(U))
V <- matrix(0, p, 10); V[] <- rnorm(length(V))
X <- tcrossprod(U, V) + 5 * matrix(rnorm(n * p), n, p)

K <- 15
system.time(svd2 <- PRIMME::svds(X, K, isreal = TRUE))

ind <- -(1:100)
X2 <- X[, ind]
u0 <- svd2$u
v0 <- svd2$v[ind, ]
system.time(
  svd2.1 <- PRIMME::svds(X2, K, isreal = TRUE))
system.time(
  svd2.2 <- PRIMME::svds(X2, K, isreal = TRUE, u0 = u0, v0 = v0))
system.time(
  svd2.3 <- PRIMME::svds(X2, K, isreal = TRUE, u0 = u0, v0 = svd(v0)$u))
system.time(
  svd2.4 <- PRIMME::svds(X2, K, isreal = TRUE, u0 = u0,
                         v0 = crossprod(X2, sweep(u0, 2, svd2$d, '/'))))
system.time(
  svd2.5 <- PRIMME::svds(X2, K, isreal = TRUE, u0 = svd2.1$u, v0 = svd2.1$v,
                         printLevel = 5))

plot(svd2.5$u, svd2.1$u, pch = 20)
all.equal(svd2.5$v, svd2.1$v)

system.time(PRIMME::svds(X2, K, isreal = TRUE))
system.time(
  svd2.21 <- PRIMME::svds(X2, K, isreal = TRUE, tol = 1e-2))
system.time(
  svd2.22 <- PRIMME::svds(X2, K, isreal = TRUE,
                          u0 = svd2.21$u, v0 = svd2.21$v))
plot(svd2.22$u, svd2.21$u, pch = 20)
system.time(
  svd2.2 <- PRIMME::svds(X2, K, isreal = TRUE, tol = 1e-2,
                         u0 = X2 %*% V[ind, ], v0 = V[ind, ]))

svd2 <- NULL
stats <- sapply(1:13, function(k) {
  svd2 <<- PRIMME::svds(X2, k, isreal = TRUE, u0 = svd2$u, v0 = svd2$v)
  svd2$stats
})
numMatvecs <- unlist(stats[1, ])
elapsedTime <- unlist(stats[3, ])
plot(numMatvecs)
plot(numMatvecs[1:11])
plot(elapsedTime[1:11])

plot(time, pch = 20)
plot(cumsum(time), pch = 20)

lapply(2:15, function(k) {
  str(cpt <- changepoint::cpt.mean(head(numMatvecs, k),
                                   penalty="Manual",pen.value=0.8,method="AMOC",test.stat="CUSUM",
                                   # penalty="SIC",method="AMOC",
                                   # Q = 1, method = "BinSeg",
                                   # penalty = "Manual",
                                   # pen.value = 0.0001,
                                   class = FALSE))
})


str(cpt <- changepoint::cpt.mean(head(numMatvecs, 4),
                                 # penalty="Manual",pen.value=0.0001,
                                 method="AMOC",#test.stat="CUSUM",
                                 # penalty="SIC",method="AMOC",
                                 # Q = 1, method = "BinSeg",
                                 # penalty = "Manual",
                                 # pen.value = 0.0001,
                                 class = FALSE))
print(cpt[2], digits = 20)

sapply(3:15, function(k) {
  changepoint::cpt.mean(head(numMatvecs, k), method = "AMOC", class = FALSE)
})
sapply(4:15, function(k) {
  changepoint::cpt.mean(head(numMatvecs, k), method = "AMOC", class = FALSE)
})
