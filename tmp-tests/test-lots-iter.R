library(bigsnpr)
obj.bed <- bed("../POPRES_data/POPRES_allchr.bed"); nPC <- 10
# obj.bed <- bed("../POPRES_data/POPRES_allchr_QC_norel.bed"); nPC <- 10
obj.svd <- bed_autoSVD2(obj.bed, k = nPC, ncores = nb_cores())
plot(obj.svd)
attr(obj.svd, "lrldr")
plot(obj.svd, type = "loadings", loadings = 3:8, coef = 0.5)
plot(obj.svd, type = "scores", scores = 3:10, coef = 0.5)

clust <- dbscan::hdbscan(obj.svd$u, 3)$cluster

plot_grid(plotlist = lapply(1:(nPC-1), function(k) {
  qplot(obj.svd$u[, k], obj.svd$u[, k+1], color = as.factor(clust)) +
    scale_colour_viridis_d() +
    coord_equal() +
    theme(legend.position = "none") +
    theme_bigstatsr(0.5)
}))

ldof <- DDoutlier::LDOF(obj.svd$u, k = 20)
lof <- LOF(obj.svd$u, seq_k = 4:15, log = FALSE)
plot(lof, ldof)
hist(sqrt(ldof), "FD"); abline(v = print(tukey_mc_up(sqrt(ldof))), col = "red")
hist(log(ldof), "FD"); abline(v = print(tukey_mc_up(log(ldof))), col = "red")

plot_grid(plotlist = lapply(1:(nPC-1), function(k) {
  qplot(obj.svd$u[, k], obj.svd$u[, k+1], color = ldof) +
    scale_colour_viridis_c() +
    coord_equal() +
    theme(legend.position = "none") +
    theme_bigstatsr(0.5)
}))

plot_grid(plotlist = lapply(1:(nPC-1), function(k) {
  qplot(obj.svd$u[, k], obj.svd$u[, k+1], color = sqrt(ldof) > tukey_mc_up(sqrt(ldof))) +
    scale_colour_viridis_d() +
    coord_equal() +
    theme(legend.position = "none") +
    theme_bigstatsr(0.5)
}))
