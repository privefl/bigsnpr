library(bigsnpr)
library(bigutilsr)
obj.bed <- bed("../paper2-PRS/backingfiles/celiacQC.bed"); nPC <- 15
# obj.bed <- bed("../POPRES_data/POPRES_allchr.bed"); nPC <- 10
# obj.bed <- bed("tmp-data/1000G_phase3_common_norel.bed"); nPC <- 20
# obj.bed <- bed(system.file("extdata", "example.bed", package = "bigsnpr")); nPC <- 10

stats <- bigsnpr:::bed_stats(obj.bed, rows_along(obj.bed), cols_along(obj.bed))
af <- stats$sum / (2 * stats$nb_nona_col)
hist(maf <- pmin(af, 1 - af))


ind.keep <- bed_clumping(obj.bed, ncores = nb_cores(), exclude = which(maf < 0.01))

obj.svd <- bed_randomSVD(obj.bed, k = nPC, ind.col = ind.keep, ncores = nb_cores())


lof <- LOF(obj.svd$u, seq_k = 3:5, log = FALSE)
hist(lof, "FD")
hist(log(lof), "FD")
hist(sqrt(lof), "FD")
hist(-(1 / lof), "FD")

library(ggplot2)
all_lof <- sapply(seq_len(nPC - 1), function(k) {
  LOF(obj.svd$u[, k + 0:1], seq_k = 3:5, log = FALSE)
})
apply(all_lof, 2, range)
hist(all_lof, "FD")
lof2 <- log(apply(all_lof, 1, max))
hist(lof2, "FD"); abline(v = print(tukey_mc_up(lof2)), col = "red")
hist(sqrt(lof2), "FD"); abline(v = print(tukey_mc_up(sqrt(lof2))), col = "red")
hist(log(lof2), "FD"); abline(v = print(tukey_mc_up(log(lof2))), col = "red")
plot(lof, lof2)

plot_grid(plotlist = lapply(1:(nPC-1), function(k) {
  qplot(obj.svd$u[, k], obj.svd$u[, k+1], color = lof) +
    scale_colour_viridis_c() +
    coord_equal() +
    theme(legend.position = "none") +
    theme_bigstatsr(0.5)
}))

plot_grid(plotlist = lapply(1:(nPC-1), function(k) {
  qplot(obj.svd$u[, k], obj.svd$u[, k+1], color = sqrt(lof2)) +
    scale_colour_viridis_c() +
    coord_equal() +
    theme(legend.position = "none") +
    theme_bigstatsr(0.5)
}))

hist(lof2, nclass.scottRob, ylim = c(0, 30)); abline(v = print(q <- tukey_mc_up(lof2, alpha = 0.05)), col = "red")
sum(out <- (lof2 > q))

plot_grid(plotlist = lapply(1:(nPC-1), function(k) {
  qplot(obj.svd$u[, k], obj.svd$u[, k+1], color = out) +
    scale_colour_viridis_d() +
    coord_equal() +
    theme(legend.position = "none") +
    theme_bigstatsr(0.5)
}))

Rcpp::sourceCpp('tmp-tests/kNN.cpp')
U <- obj.svd$u
u1 <- U[, 1]
ord <- order(u1)
u2 <- U[, 2]

knn0 <- bigutilsr::knn(U[, 1:2], k = 11)$nn.dists[, -1]
knn <- dist_kNN(u1[ord], U[ord, 2], nfirst = 10)
plot(knn, knn0[ord, ])
microbenchmark::microbenchmark(
  NABOR = bigutilsr::knn(U[, 1:2], k = 11)$nn.dists[, -1],
  ME = {
    u1 <- U[, 1]
    ord <- order(u1)
    u2 <- U[, 2]
    dist_kNN(u1[ord], U[ord, 2], nfirst = 10)
  },
  times = 5
)
# Unit: microseconds
#  expr     min      lq     mean  median      uq     max neval
# NABOR 610.371 620.173 692.1248 644.453 774.784 810.843     5
#    ME 580.355 593.685 602.6050 606.893 608.618 623.474     5
# Unit: milliseconds
#  expr       min        lq      mean   median        uq       max neval
# NABOR  22.85622  23.10987  24.23842  23.6577  25.62727  25.94106     5
#    ME 104.96057 105.54911 108.88959 106.2115 110.43133 117.29546     5

round(bigutilsr::knn(obj.svd$u, k = 11)$nn.dists[order(lof2, decreasing = TRUE)[1:5], -1], 4)

dist <- bigutilsr::knn(obj.svd$u, k = 11)$nn.dists[, -1]
plot(dist[, 2:3+5])

dist2 <- dist %*% (1/(1:10))
plot_grid(plotlist = lapply(1:(nPC-1), function(k) {
  qplot(obj.svd$u[, k], obj.svd$u[, k+1], color = dist2) +
    scale_colour_viridis_c() +
    coord_equal() +
    theme(legend.position = "none") +
    theme_bigstatsr(0.5)
}))
hist(dist2)
