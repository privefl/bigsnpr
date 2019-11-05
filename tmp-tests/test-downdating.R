library(bigsnpr)
library(bigutilsr)
# obj.bed <- bed("../paper2-PRS/backingfiles/celiacQC.bed"); nPC <- 20
obj.bed <- bed("../POPRES_data/POPRES_allchr.bed"); nPC <- 15
# obj.bed <- bed("tmp-data/1000G_phase3_common_norel.bed"); nPC <- 20

stats <- bigsnpr:::bed_stats(obj.bed, rows_along(obj.bed), cols_along(obj.bed))
af <- stats$sum / (2 * stats$nb_nona_col)
hist(maf <- pmin(af, 1 - af))


ind.keep <- bed_clumping(obj.bed, ncores = nb_cores(), exclude = which(maf < 0.01))

obj.svd <- bed_randomSVD(obj.bed, k = nPC, ind.col = ind.keep, ncores = nb_cores())
stopifnot(all.equal(2 * af[ind.keep], obj.svd$center))

# plot(obj.svd, type = "loadings", loadings = 1:12, coeff = 0.5)

library(bigutilsr)
S.col <- log(rowSums(apply(obj.svd$v, 2, function(x) (x - median(x))^2)))
# roll mean to get only consecutive outliers
hist(S2.col <- rollmean(S.col, size = 50),
     breaks = nclass.scottRob, ylim = c(0, 100))
abline(v = tukey_mc_up(S2.col), col = "red")

length(ind <- -which(S2.col > tukey_mc_up(S2.col)))

A <- bigsnpr:::read_bed_scaled(
  obj.bed, rows_along(obj.bed), ind.keep[-ind],
  center = obj.svd$center[-ind],
  scale = -obj.svd$scale[-ind]
)
UT_A <- crossprod(obj.svd$u, A)
M_A <- A - obj.svd$u %*% UT_A
svd_A <- svd(M_A)
U_A <- svd_A$u

BT_V <- obj.svd$v[-ind, ]
eig_B <- eigen(diag(nrow(BT_V)) - tcrossprod(BT_V), symmetric = TRUE)
VD_B <- sweep(eig_B$vectors, 2, sqrt(eig_B$values), '*')
VDI_B <- sweep(VD_B, 2, eig_B$values, '/')

Q1 <- UT_A %*% BT_V + diag(obj.svd$d)
VD_A <- sweep(svd_A$v, 2, svd_A$d, '*')
Q2 <- crossprod(VD_A, BT_V)
# VD_B <- sweep(svd_B$v, 2, svd_B$d, '*')
Q3 <- UT_A %*% VD_B
Q4 <- crossprod(VD_A, VD_B)
Q <- cbind(rbind(Q1, Q2), rbind(Q3, Q4))
svd_Q <- svd(Q, nu = ncol(obj.svd$u), nv = ncol(obj.svd$v))
U2 <- cbind(obj.svd$u, U_A) %*% svd_Q$u
BT_V2 <- obj.svd$v[ind, ]
V2 <- cbind(BT_V2, -BT_V2 %*% crossprod(BT_V, VDI_B)) %*% svd_Q$v

round(100 * cor(U2, obj.svd$u)^2, 1)
round(100 * cor(V2, obj.svd$v[ind, ])^2, 1)
sum(cor(U2, obj.svd$u)^2 > 0.9)
sum(cor(V2, obj.svd$v[ind, ])^2 > 0.9)
# plot(obj.svd$u)
# plot(U2)
round(corr <- 100 * cor(U2, obj.svd$u) * cor(V2, obj.svd$v[ind, ]), 1)
sum(corr > 95) # 7 > 95 & 9 > 90
which(corr > 95, arr.ind = TRUE)

system.time(
  obj.svd2 <- bed_randomSVD(obj.bed, k = nPC, ind.col = ind.keep[ind],
                            ncores = nb_cores())
) # 700 sec (sequential)
round(100 * cor(U2, obj.svd2$u), 1)
round(100 * cor(V2, obj.svd2$v), 1)
round(100 * crossprod(U2), 3)
round(100 * crossprod(V2), 3)
all.equal(V2[, 7], rep(0, nrow(V2)))
plot(svd_Q$d, log = "xy")
plot(V2[, 7], pch = 20)
plot(V2[, 8], pch = 20)
plot(U2[, 7:8])

# system.time(
#   obj.svd2 <- bed_randomSVD(obj.bed, k = nPC, ind.col = ind.keep[ind],
#                             v0 = V2[, seq_len(sum(corr > 95)), drop = FALSE])
# ) # 770 with first 7 PC loadings // 663

plot(obj.svd2$u[, 5:6 + 4])

S.col <- log(rowSums(apply(V2, 2, function(x) (x - median(x))^2)))
# roll mean to get only consecutive outliers
hist(S2.col <- rollmean(S.col, size = 50),
     breaks = nclass.scottRob, ylim = c(0, 100))
abline(v = tukey_mc_up(S2.col), col = "red")

S.col <- log(rowSums(apply(obj.svd2$v, 2, function(x) (x - median(x))^2)))
# roll mean to get only consecutive outliers
hist(S2.col <- rollmean(S.col, size = 50),
     breaks = nclass.scottRob, ylim = c(0, 100))
abline(v = tukey_mc_up(S2.col), col = "red")
