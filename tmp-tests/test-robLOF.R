library(bigsnpr)
library(bigutilsr)
# obj.bed <- bed("../paper2-PRS/backingfiles/celiacQC.bed"); nPC <- 8
# obj.bed <- bed("../POPRES_data/POPRES_allchr.bed"); nPC <- 8
obj.bed <- bed("tmp-data/1000G_phase3_common_norel.bed"); nPC <- 20

stats <- bigsnpr:::bed_stats(obj.bed, rows_along(obj.bed), cols_along(obj.bed))
af <- stats$sum / (2 * stats$nb_nona_col)
hist(maf <- pmin(af, 1 - af))


ind.keep <- bed_clumping(obj.bed, ncores = nb_cores(), exclude = which(maf < 0.01))

obj.svd <- bed_randomSVD(obj.bed, k = nPC, ind.col = ind.keep, ncores = nb_cores())

all_lof <- sapply(c(4, 10, 30), function(kNN) {
  bigutilsr::LOF(obj.svd$u, seq_k = kNN, log = FALSE)
})
new_lof <- bigutilsr::LOF(obj.svd$u, seq_k = 4:8, log = FALSE,
                          combine = max)#function(x) log(mean(exp(x))))
hist(new_llof <- log(new_lof), "FD")
new_lof2 <- bigsnpr:::robLOF(obj.svd$u, seq_kNN = 4:8, log = FALSE,
                          combine = max,#function(x) log(mean(exp(x))),
                          ncores = nb_cores())
hist(new_llof2 <- log(new_lof2), "FD")
abline(v = tukey_mc_up(new_llof2), col = "red")
plot(new_lof, new_lof2); abline(0, 1, col = "red")

new_lof3 <- bigutilsr::LOF(obj.svd$u, seq_k = 4:8, log = FALSE, robMaha = TRUE,
                          combine = max)#function(x) log(mean(exp(x))))
hist(new_llof3 <- log(new_lof3), "FD")
plot(new_lof, new_lof3); abline(0, 1, col = "red")

U <- obj.svd$u
maha <- covRob(U, estim = "pairwiseGK", distance = FALSE, corr = FALSE)
eigs <- eigen(unname(maha$cov), symmetric = TRUE)
U <- U %*% sweep(eigs$vectors, 2, sqrt(eigs$values),
                 "/")
plot(U[, 3:4 + 14])

plot_grid(plotlist = lapply(1:10, function(k) {
  plot(obj.svd, type = "scores", scores = 2 * k - (1:0)) +
    aes(color = new_llof3) +
    scale_colour_viridis_c() + theme(legend.position = "none")
}))
