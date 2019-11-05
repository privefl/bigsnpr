library(bigsnpr)
library(bigutilsr)
# obj.bed <- bed("../paper2-PRS/backingfiles/celiacQC.bed"); nPC <- 8
obj.bed <- bed("../POPRES_data/POPRES_allchr.bed"); nPC <- 16
# obj.bed <- bed("tmp-data/1000G_phase3_common_norel.bed"); nPC <- 20

stats <- bigsnpr:::bed_stats(obj.bed, rows_along(obj.bed), cols_along(obj.bed))
af <- stats$sum / (2 * stats$nb_nona_col)
hist(maf <- pmin(af, 1 - af))


ind.keep <- bed_clumping(obj.bed, ncores = nb_cores(), exclude = which(maf < 0.01))

obj.svd <- bed_randomSVD(obj.bed, k = nPC, ind.col = ind.keep, ncores = nb_cores())


knn <- nabor::knn(predict(obj.svd), k = 6)
dists <- knn$nn.dists[, -1, drop = FALSE]
dist_mean <- rowMeans(dists^2)
hist(s <- log(dist_mean), nclass.scottRob)
abline(v = print(tukey_mc_up(s)), col = "red")
abline(v = print(q <- hist_out(s, nboot = 100)$lim[2]), col = "blue")

library(ggplot2)
plot_grid(plotlist = lapply(1:(nPC/2), function(k) {
  K <- 2 * k
  qplot(obj.svd$u[, K - 1], obj.svd$u[, K], color = s2) +
    scale_colour_viridis_c() +
    theme_bigstatsr() +
    theme(legend.position = "none") +
    labs(x = paste("PC", K - 1), y = paste("PC", K))
}))
plot_grid(plotlist = lapply(1:(nPC/2), function(k) {
  K <- 2 * k
  qplot(obj.svd$u[, K - 1], obj.svd$u[, K], color = (s > q)) +
    scale_colour_viridis_d() +
    theme_bigstatsr() +
    theme(legend.position = "none") +
    labs(x = paste("PC", K - 1), y = paste("PC", K))
}))


obj.svd2 <- bed_randomSVD(obj.bed, k = nPC, ind.row = which(!(s > q)),
                          ind.col = ind.keep, ncores = nb_cores())

knn2 <- nabor::knn(predict(obj.svd2), k = 6)
dists2 <- knn2$nn.dists[, -1, drop = FALSE]
dist_mean2 <- rowMeans(dists2^2)
hist(s2 <- log(dist_mean2), nclass.scottRob)
abline(v = print(tukey_mc_up(s2)), col = "red")
abline(v = print(q2 <- hist_out(s2, nboot = 100)$lim[2]), col = "blue")

library(ggplot2)
plot_grid(plotlist = lapply(1:(nPC/2), function(k) {
  K <- 2 * k
  qplot(obj.svd2$u[, K - 1], obj.svd2$u[, K], color = s2) +
    scale_colour_viridis_c() +
    theme_bigstatsr() +
    theme(legend.position = "none") +
    labs(x = paste("PC", K - 1), y = paste("PC", K))
}))
plot_grid(plotlist = lapply(1:(nPC/2), function(k) {
  K <- 2 * k
  qplot(obj.svd2$u[, K - 1], obj.svd2$u[, K], color = (s2 > q2)) +
    scale_colour_viridis_d() +
    theme_bigstatsr() +
    theme(legend.position = "none") +
    labs(x = paste("PC", K - 1), y = paste("PC", K))
}))
