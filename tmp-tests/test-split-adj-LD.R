library(bigsnpr)
library(dplyr)
library(ggplot2)

obj.bed <- bed(download_1000G(dir = "../datasets"))
fam2 <- bigreadr::fread2(paste0(obj.bed$prefix, ".fam2"))
THR_R2 <- 0.02
CHR <- 22

# plink2 <- download_plink2("../datasets")
# rel <- snp_plinkKINGQC(plink2, obj.bed$bedfile, thr.king = 0.044,
#                        make.bed = FALSE, ncores = 4)
# dput(which(obj.bed$fam$sample.ID %in% rel$IID2))
ind_rel <- c(19L, 80L, 212L, 238L, 240L, 310L, 619L, 900L, 902L, 1005L, 1037L,
             1176L, 1263L, 1269L, 1271L, 1291L, 1308L, 1309L, 1314L, 1493L, 1619L,
             1682L, 1978L, 1986L, 2084L, 2095L, 2104L, 2108L, 2171L, 2173L, 2405L)

# prepare for snp_splitld
ind.row <-  setdiff(which(fam2$`Super Population` == "EAS"), ind_rel)
ind.col <- which(obj.bed$map$chromosome == CHR)
obj.bigsnp <- snp_attach(
  snp_readBed2(obj.bed$bedfile, tempfile(), ind.row, ind.col, ncores = 4))
G <- obj.bigsnp$genotypes
ind.col2 <- which(snp_MAF(G, ncores = 4) > 0.05)
map <- obj.bigsnp$map[ind.col2, ]
POS <- map$physical.pos
POS2 <- snp_asGeneticPos(map$chromosome, POS, dir = "../datasets")

corr <- snp_cor(G, ind.col = ind.col2, infos.pos = POS2, size = 3 / 1000,
                thr_r2 = THR_R2, ncores = 4)
dim(corr)
corr[1:5, 1:5]

corr0 <- snp_cor_with_cov(G, ind.col = ind.col2, covar.row = G[, 0],
                          infos.pos = POS2, size = 3 / 1000,
                          thr_r2 = THR_R2, ncores = 4)
all.equal(corr0, corr)

# snp_splitld
res <- bigsnpr::snp_ldsplit(
  corr, thr_r2 = THR_R2, min_size = 200, max_size = 3000, max_K = 20)

qplot(data = res, n_block, cost) +
  theme_bw(12) +
  scale_x_continuous(breaks = 0:100 * 10, minor_breaks = 0:500 * 2) +
  scale_y_log10() +
  theme(legend.position = "top") +
  labs(y = "Sum of squared correlations outside blocks (log-scale)", x = "Number of blocks")

## With covariates
obj.svd <- snp_autoSVD(G, obj.bigsnp$map$chromosome, ncores = 4)
plot(obj.svd)
plot(obj.svd, type = "scores", scores = 1:8)

corr2 <- snp_cor_with_cov(G, ind.col = ind.col2, covar.row = obj.svd$u[, 1:2],
                          infos.pos = POS2, size = 3 / 1000,
                          thr_r2 = THR_R2, ncores = 4)
corr2
corr2[1:5, 1:5]
ind <- Matrix::which(corr != 0 & corr2 != 0)
ind2 <- sort(sample(ind, 50e3))
plot(corr[ind2], corr2[ind2]); abline(0, 1, col = "red", lwd = 2)


res2 <- bigsnpr::snp_ldsplit(
  corr2, thr_r2 = THR_R2, min_size = 200, max_size = 3000, max_K = 20)

qplot(data = res2, n_block, cost) +
  geom_point(data = res, color = "red") +
  theme_bw(12) +
  scale_x_continuous(breaks = 0:100 * 10, minor_breaks = 0:500 * 2) +
  scale_y_log10() +
  theme(legend.position = "top") +
  labs(y = "Sum of squared correlations outside blocks (log-scale)", x = "Number of blocks")


obj.svd2 <- bed_autoSVD(obj.bed, ind.row, ncores = 4)
plot(obj.svd2)
plot(obj.svd2, type = "scores", scores = 1:8)

corr3 <- snp_cor_with_cov(G, ind.col = ind.col2, covar.row = obj.svd2$u[, 1:5],
                          infos.pos = POS2, size = 3 / 1000,
                          thr_r2 = THR_R2, ncores = 4)

res3 <- bigsnpr::snp_ldsplit(
  corr3, thr_r2 = THR_R2, min_size = 200, max_size = 3000, max_K = 20)

qplot(data = res2, n_block, cost) +
  geom_point(data = res, color = "red") +
  geom_point(data = res3, color = "blue") +
  theme_bw(12) +
  scale_x_continuous(breaks = 0:100 * 10, minor_breaks = 0:500 * 2) +
  scale_y_log10() +
  theme(legend.position = "top") +
  labs(y = "Sum of squared correlations outside blocks (log-scale)", x = "Number of blocks")

ind <- Matrix::which(corr != 0 | corr3 != 0)
ind2 <- sort(sample(ind, 50e3))
plot(corr[ind2], corr3[ind2]); abline(0, 1, col = "red", lwd = 2)
mean(corr[ind2]^2 > corr3[ind2]^2)  # 57.2%
mean(corr[ind2]^2 - corr3[ind2]^2)  # 0.002
