library(bigsnpr)
# obj.bed <- bed("../paper2-PRS/backingfiles/celiacQC.bed"); nPC <- 20
obj.bed <- bed("../POPRES_data/POPRES_allchr.bed"); nPC <- 10
# obj.bed <- bed("tmp-data/1000G_phase3_common_hapmap_norel.bed"); nPC <- 20

obj.svd <- bed_autoSVD2(obj.bed, k = nPC, ncores = nb_cores())

plot(obj.svd)
attr(obj.svd, "lrldr")
plot(obj.svd, type = "loadings", loadings = 3:8, coef = 0.5)
plot(obj.svd, type = "scores", scores = 3:10, coef = 0.5)

gwas <- bed_pcadapt(obj.bed, obj.svd$u, ind.row = attr(obj.svd, "subset.row"),
                    ind.col = attr(obj.svd, "subset.col"), ncores = nb_cores())
hist(sqrt(gwas$score), "FD")
plot(-predict(gwas), pch = 20)

V <- obj.svd$v
maha <- covRob(V, estim = "pairwiseGK", distance = FALSE, corr = FALSE)
eigs <- eigen(unname(maha$cov), symmetric = TRUE)
V <- V %*% sweep(eigs$vectors, 2, sqrt(eigs$values), "/")
score2 <- rowSums(V^2)
score3 <- stats::mahalanobis(obj.svd$v, center = maha$center, cov = maha$cov)
score4 <- bigutilsr::covRob(obj.svd$v, corr = FALSE, estim = "pairwiseGK")$dist
plot(gwas$score, score2)
plot(score3, score2)
plot(gwas$score, score3)
plot(score3, score4)
