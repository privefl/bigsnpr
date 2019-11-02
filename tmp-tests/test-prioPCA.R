library(bigsnpr)
library(bigutilsr)
obj.bed <- bed("../paper2-PRS/backingfiles/celiacQC.bed"); nPC <- 10
# obj.bed <- bed("../POPRES_data/POPRES_allchr.bed"); nPC <- 10
# obj.bed <- bed(download_1000G("~/Desktop/paper4-bedpca/tmp-data")); nPC <- 20
# obj.bed <- bed(system.file("extdata", "example.bed", package = "bigsnpr")); nPC <- 10

stats <- bigsnpr:::bed_stats(obj.bed, rows_along(obj.bed), cols_along(obj.bed))
af <- stats$sum / (2 * stats$nb_nona_col)
hist(maf <- pmin(af, 1 - af))
is_rm_for_good <- (maf < 0.01)
keep <- rep(FALSE, ncol(obj.bed))
obj.pcadapt <- list()

ind.keep <- bed_clumping(obj.bed, S = obj.pcadapt$score, ncores = nb_cores(),
                         exclude = which(is_rm_for_good))
keep[ind.keep] <- TRUE

obj.svd <- bed_randomSVD(obj.bed, k = nPC, ind.col = which(keep), ncores = nb_cores())
plot(obj.svd)
plot(obj.svd, type = "scores", scores = 1:8, coeff = 0.5)

obj.pcadapt <- bed_pcadapt(obj.bed, obj.svd$u, ind.col = which(keep), ncores = nb_cores())
plot(obj.pcadapt)

hist(S2 <- rollmean(S <- sqrt(obj.pcadapt$score), 50), "FD")
abline(v = print(q <- tukey_mc_up(S2)), col = "red")
investigate <- which(keep)[S2 > q]

S3 <- rep(NA, ncol(obj.bed))
S3[keep] <- S

ind.keep2 <- bed_clumping(obj.bed, S = S3,
                          ncores = nb_cores(),
                          exclude = !cols_along(obj.bed) %in% investigate)
keep[investigate] <- FALSE
keep[ind.keep2] <- TRUE
sum(keep)
obj.svd2 <- bed_randomSVD(obj.bed, k = nPC, ind.col = ind.keep2, ncores = nb_cores())
plot(obj.svd)
plot(obj.svd2)
plot(obj.svd2, type = "scores", scores = 1:8, coeff = 0.5)
round(cor(obj.svd2$u, obj.svd$u), 2)
plot(obj.svd, type = "scores", scores = c(4, 6))
plot(obj.svd2, type = "scores", scores = c(4, 6))


ind.keep3 <- bed_clumping(obj.bed, S = obj.pcadapt$score, thr.r2 = 0.1,
                          ncores = nb_cores(), exclude = which(S < q))
ind.keep4 <- c(setdiff(ind.keep2, which(S > q)), ind.keep3)
obj.svd3 <- bed_randomSVD(obj.bed, k = nPC, ind.col = ind.keep4, ncores = nb_cores())
plot(obj.svd3)
plot(obj.svd3, type = "scores", scores = 1:8, coeff = 0.5)

obj.pcadapt2 <- bed_pcadapt(obj.bed, obj.svd3$u, ncores = nb_cores())
plot(obj.pcadapt2)
hist(S <- rollmean(sqrt(obj.pcadapt2$score), 50), "FD")
abline(v = print(q <- tukey_mc_up(S)), col = "red")
ind.keep5 <- bed_clumping(obj.bed, S = obj.pcadapt2$score, thr.r2 = 0.01,
                          ncores = nb_cores(), exclude = which(S < q))
ind.keep6 <- c(setdiff(ind.keep3, which(S > q)), ind.keep5)
