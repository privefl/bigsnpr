k1 <- k2 <- 10
knn <- nabor::knn(obj.svd$u, k = max(k1, k2) + 1)
ids <- knn$nn.idx[, -1, drop = FALSE]
dists <- knn$nn.dists[, -1, drop = FALSE]


mean_d <- rowMeans(dists[, 1:k1]^2)
lof.num2 <- sapply(seq_len(nrow(ids)), function(i) {
  mean(mean_d[ids[i, 1:k2]])
})

library(bigutilsr)
# hist(lof <- sqrt(LOF(obj.svd$u, seq_k = k1, log = FALSE)), "FD")
hist(sqrt(mean_d), "FD")
abline(v = print(tukey_mc_up(sqrt(mean_d))), col = "red")
hist(lof2 <- log(mean_d / lof.num2), "FD")
abline(v = print(tukey_mc_up(lof2)), col = "red")
hist(lof3 <- mean_d / sqrt(lof.num2), "FD")
abline(v = print(tukey_mc_up(lof3)), col = "red")
hist(lof4 <- log(lof3), "FD")
abline(v = print(tukey_mc_up(lof4)), col = "red")
# hist(lof5 <- log(mean_d * (lof.num2)^(1/3)), "FD")
# abline(v = print(tukey_mc_up(lof5)), col = "red")

PLOF <- function (dataset, k = 5) {
  n <- nrow(dataset)
  dist.obj <- dbscan::kNN(dataset, k)
  nnSD <- apply(dist.obj$dist, 1, function(x) {
    sqrt((sum(x^2)/k))
  })
  plof <- NULL
  for (i in 1:n) {
    plof[i] <- (nnSD[i]/mean(nnSD[dist.obj$id[i, ]])) - 1
  }
  plof
}

plot(lof2, PLOF(obj.svd$u, k = k1))

PLOF2 <- function(U, kNN, ncores) {
  all_knn <- big_parallelize(U, function(X, ind, k) {
    bigutilsr::knn(X, X[ind, , drop = FALSE], k = k)
  }, ind = rows_along(U), k = kNN + 1, ncores = ncores)
  dists <- do.call("rbind", lapply(all_knn, function(x) {
    x$nn.dists[, -1, drop = FALSE]
  }))
  plof.num <- sqrt(rowMeans(dists^2))
  ids <- do.call("rbind", lapply(all_knn, function(x) {
    x$nn.idx[, -1, drop = FALSE]
  }))
  plof.deno <- sapply(rows_along(ids), function(i) {
    mean(plof.num[ids[i, 1:k2]])
  })
  plof.num / plof.deno - 1
}
plot(plof <- PLOF2(obj.svd$u, k = k1), PLOF(obj.svd$u, k = k1))
hist(log1p(plof), "FD")
abline(v = print(tukey_mc_up(log1p(plof))), col = "red")

plof.num <- sqrt(rowMeans(dists^2))
plof.deno <- sapply(rows_along(ids), function(i) {
  mean(plof.num[ids[i, 1:k2]])
})
hist(plof <- plof.num / sqrt(plof.deno), "FD")

hist(log(plof))
hist(plof2 <- -1 / (plof), "FD")
abline(v = print(q <- tukey_mc_up(plof2, alpha = 0.05)), col = "red")
sum(plof2 > q)

hist(lof <- bigutilsr::LOF(predict(obj.svd), seq_k = 4:6, log = FALSE))
hist(lof2 <- -1 / lof, "FD")
abline(v = print(q <- tukey_mc_up(lof2, alpha = 0.05)), col = "red")
sum(lof2 > q)


library(ggplot2)
plot_grid(plotlist = lapply(1:(nPC-1), function(k) {
  qplot(obj.svd$u[, k], obj.svd$u[, k+1], color = plof) +
    scale_colour_viridis_c() +
    coord_equal() +
    theme_bigstatsr(0.5) +
    theme(legend.position = "none")
}))

plot_grid(plotlist = lapply(1:(nPC-1), function(k) {
  qplot(obj.svd$u[, k], obj.svd$u[, k+1], color = (lof2)) +
    scale_colour_viridis_c() +
    coord_equal() +
    theme_bigstatsr(0.5) +
    theme(legend.position = "none")
}))

plot_grid(plotlist = lapply(1:(nPC-1), function(k) {
  qplot(obj.svd$u[, k], obj.svd$u[, k+1], color = lof4) +
    scale_colour_viridis_c() +
    coord_equal() +
    theme_bigstatsr(0.5) +
    theme(legend.position = "none")
}))

all_thr <- replicate(100, tukey_mc_up(sample(lof3, replace = TRUE)))
hist(lof3, nclass.scottRob)
abline(v = print(median(all_thr)), col = "blue")
abline(v = print(quantile(all_thr, 0.9)), col = "red")

abline(v = print(q <- hist_out(lof3, nboot = 100)$lim[2]), col = "green")

plot_grid(plotlist = lapply(1:(nPC-1), function(k) {
  qplot(obj.svd$u[, k], obj.svd$u[, k+1], color = lof3 > q) +
    theme_bigstatsr(0.5) +
    scale_colour_viridis_d() +
    coord_equal() +
    theme(legend.position = "none")
}))

hist(lof3, "FD")
mix1 <- mixtools::normalmixEM(hist_out(lof4)$x, lambda = 0.5,
                              mu = quantile(lof3, c(1, 3) / 4))

plot(mix1, whichplots = 2)
alpha <- 0.05 / length(lof3)
q <- mix1$lambda[1] * qnorm(alpha, mean = mix1$mu[1], sd = mix1$sigma[1], lower.tail = FALSE) +
  mix1$lambda[2] * qnorm(alpha, mean = mix1$mu[2], sd = mix1$sigma[2], lower.tail = FALSE)
hist(lof4, "FD")
abline(v = print(q), col = "red")
abline(v = print(tukey_mc_up(lof4)), col = "blue")

plot(MASS::fitdistr(lof4, "normal"))
