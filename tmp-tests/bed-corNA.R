library(bigsnpr)
library(bigreadr)
library(dplyr)
library(ggplot2)

rel <- fread2("data/ukb25589_rel_s488346.dat")
fam <- fread2("data/ukbb_bed/ukbb_488282.fam")
ind.row <- which(!fam$V2 %in% rel$ID2)
obj.bed <- bed("data/ukbb_bed/ukbb_488282.bed")

obj.svd <- readRDS("tmp-results/SVD_UKBB4.rds")
ind.keep <- attr(obj.svd, "subset")

U <- obj.svd$u
system.time(test <- bigsnpr:::bed_corNA(obj.bed, ind.row, ind.keep, U))
# 4 minutes
str(test)

nb_na <- test$nb_na
sumU <- colSums(U)
sumUU <- colSums(U^2)
n <- length(ind.row)
r <- (n * test$prod - outer(sumU, nb_na)) /
  outer(sqrt(n * sumUU - sumU^2), sqrt(nb_na * (n - nb_na)))
hist(r)
pval <- 2 * pt(abs(r) * sqrt((n - 2) / (1 - r^2)), df = n - 2, lower.tail = FALSE)
hist(pval)
hist(-log10(pval), xlim = c(0, 10), "FD")
hist(-log10(pval[16, ]), xlim = c(0, 10), "FD")

alpha <- 0.05 / sum(!is.na(pval))
ind.corna <- unique(which(pval < alpha, arr.ind = TRUE)[, 2])

plot(nb_na, -log10(pval[16, ]), col = scales::alpha("black", 0.3), pch = 20)
test$prod[, 1:5]
pval[, 1:5]

cor(obj.svd$v, test$nb_na)

plot(obj.svd$v[, 16], r[16, ], pch = 20, col = scales::alpha("black", 0.3))

d <- bigutilsr::covRob(cbind(obj.svd$v[, 16], r[16, ]), estim = "pairwiseGK")$dist
qplot(obj.svd$v[, 16], r[16, ], color = d) +
  scale_color_viridis_c(trans = "log") +
  theme_bigstatsr()



stats <- bigsnpr:::bed_stats(obj.bed, ind.row, ind.keep)
round(100 * cor(stats$nb_nona_row, obj.svd$u), 1)
round(100 * cor(stats$nb_nona_col, obj.svd$v), 1)

tmp <- tempfile()
system(glue::glue(
  "{download_plink('tmp-data')}",
  " --bfile {obj.bed$prefix}",
  " --test-mishap",
  " --out {tmp}"
)) # Total genotyping rate is 0.99735.
miss <- fread2(paste0(tmp, ".missing.hap"))
hist(miss$P)
ind.bad <- na.omit(match(unique(obj.bed$map$marker.ID[which(miss$P < 0.01)]),
                         obj.bed$map$marker.ID[ind.keep]))
plot(obj.svd, type = "loadings", loadings = 16) +
  aes(color = seq_along(ind.keep) %in% ind.bad)

snp_readBed2(obj.bed$bedfile, ind.row = ind.row, ind.col = ind.keep,
             backingfile = "data/UKBB_PCA")
ukbb <- snp_attach("/data/privef/UKBB_PCA.rds")
G <- ukbb$genotypes
table(G[, 1], exclude = NULL)
ind.bad <- which(predict(obj.svd)[, 16] > 50)
counts <- big_counts(G, ind.bad)
hist(pNA <- counts[4, ] / length(ind.bad))
round(100 * cor(obj.svd$u, is.na(G[, which(pNA > 0.3)])), 1)
hist(obj.svd$v[, 16]); hist(obj.svd$v[pNA > 0.3, 16], add = TRUE)
plot(pNA, obj.svd$v[, 16], pch = 20)

G.bad <- G[ind.bad, ]
ind_na <- which(is.na(G.bad), arr.ind = TRUE)

G2 <- snp_fastImputeSimple(G, method = "mean2", ncores = nb_cores())
G2[, 1]
G2[is.na(G[, 1]), 1]
obj.svd2 <- big_randomSVD(G2, snp_scaleBinom(), k = 20, ncores = nb_cores(), verbose = TRUE)
plot(obj.svd2)
round(100 * cor(obj.svd2$u, obj.svd$u), 1) # same

G3 <- snp_fastImputeSimple(G, method = "random", ncores = nb_cores())
mean(print(G3[, 1]), na.rm = TRUE)
mean(print(G3[is.na(G[, 1]), 1]))
obj.svd3 <- big_randomSVD(G3, snp_scaleBinom(), k = 20, ncores = nb_cores(), verbose = TRUE)
plot(obj.svd3)
round(100 * cor(obj.svd3$u, obj.svd$u), 1)
plot(obj.svd3, type = "scores", scores = 15:16)


snp_qc <- bigreadr::fread2("https://biobank.ctsu.ox.ac.uk/crystal/crystal/auxdata/ukb_snp_qc.txt")
names(snp_qc)
ind.col <- which(rowSums(snp_qc[match(obj.bed$map$marker.ID, snp_qc$rs_id), 10:115]) == 106)
obj.svd4 <- bed_autoSVD(obj.bed, ind.row = ind.row, ind.col = ind.col,
                        k = 20, ncores = nb_cores())
# saveRDS(obj.svd4, "tmp-results/SVD_UKBB4.rds")
plot(obj.svd4)
plot(obj.svd4, type = "loadings", loadings = 15:18, coeff = 0.7)
ind.keep2 <- attr(obj.svd4, "subset")
all(snp_qc[match(obj.bed$map$marker.ID[ind.keep2], snp_qc$rs_id), 10:115] == 1)

stats2 <- bigsnpr:::bed_stats(obj.bed, ind.row, ind.keep2)
-log10(apply(obj.svd4$v, 2, function(x) cor.test(x, stats2$nb_nona_col)$p.value))
