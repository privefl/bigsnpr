X <- matrix(rnorm(10e3 * 10), ncol = 10)
dist <- robust::covRob(X, estim = "pairwiseGK")$dist
hist(dist)
median(dist)

N <- 4e3
X2 <- rbind(X[rep(1, N), ], X)
dist2 <- robust::covRob(X2, estim = "pairwiseGK")$dist
hist(dist2)
median(dist2)
plot(dist, dist2[-(1:N)])
dist2[1:5]
dist[1]
