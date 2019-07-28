N <- 50e3  # 10e3 -> 50e3
U <- matrix(rnorm(N * 20), N)
U <- readRDS("tmp-data/svd_celiac.rds")$u

system.time(
  knn <- nabor::knn(U, k = 31, searchtype = "kd_tree_heap")
) # 3.0 -> 203
system.time(
  knn <- nabor::knn(U, k = 31, searchtype = "kd_linear_heap")
) # 3.2 -> 243
system.time(
  knn2 <- dbscan::kNN(U, k = 30)
) # 6.2 -> 311

d <- knn2
ids <- knn$nn.idx[, -1]
dists <- knn$nn.dists[, -1]
lof3 <- sapply(c(4, 10, 30), function(k) {
  print(k)

  lrd <- sapply(rows_along(dists), function(i) {
    maxs <- pmax(dists[ids[i, 1:k], k], dists[i, 1:k])
    1 / mean(maxs)
  })

  lof <- sapply(rows_along(ids), function(i) {
    mean(lrd[ids[i, 1:k]]) / lrd[i]
  })

  lof
})
hist(log(apply(lof3, 1, max)))


# all.equal(knn2$id, knn$nn.idx[, -1], check.attributes = FALSE)
all.equal(knn2$dist, knn$nn.dists[, -1], check.attributes = FALSE)
# plot(knn2$dist, knn$nn.dists[, -1])
knn2$dist[1, ]
knn$nn.dists[1, ]

system.time(
  knn3 <- RcppHNSW::hnsw_knn(U, k = 30, distance = "l2")
) # 2.5

system.time(
  knn4 <- FNN::get.knn(U, k = 30)
) # 2.8
