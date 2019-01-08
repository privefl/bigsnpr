library(bigsnpr)
prefix <- "tmp-data/indep"
plink <- download_plink()
bed2 <- snp_plinkQC(plink, prefix.in = prefix,
                    maf = 0.05, geno = 0.05, mind = 0.05, autosome.only = TRUE)
bed3 <- snp_plinkIBDQC(plink, bed2, ncores = nb_cores(), pi.hat = 0.12,
                       pruning.args = c(500, 0.1))

library(bigsnpr)
bed3 <- "tmp-data/indep_QC_norel.bed"
rds <- snp_readBed(bed3)
snp <- snp_attach(rds)
G <- snp$genotypes
CHR <- snp$map$chromosome
counts <- big_counts(G)
counts[, 1:10]

system.time(
  infos <- snp_fastImpute(G, CHR, seed = 1, ncores = nb_cores(), thr.imp = 0.01)
) # 424 for 0.01
# In the meantime, in another session, you can check progress
infos <- big_attach("tmp-data/indep_QC_norel-infos-impute.rds")
mean(!is.na(infos[1, ]))

pvals <- c(0.01, 0.005, 0.002, 0.001); colvals <- 2:5
df <- data.frame(pNA = infos[1, ], pError = infos[2, ])

# base R
plot(subset(df, pNA > 0.01), pch = 20)
curve(0.01 / x, from = 0, lwd = 2, col = "red", add = TRUE)
snp$genotypes$code256 <- CODE_IMPUTE_PRED
rds <- subset(snp, ind.col = with(df, which(pNA * pError < 0.01)))
snp2 <- snp_attach(rds)

library(bigsnpr)
snp2 <- snp_attach("tmp-data/indep_QC_norel_sub1.rds")
G <- snp2$genotypes
CHR <- snp2$map$chromosome
POS <- snp2$map$physical.pos
svd <- big_randomSVD(G, snp_scaleBinom(), k = 20, ncores = nb_cores())
plot(svd)
plot(svd, type = "scores")
plot(svd, type = "scores", scores = 9:10)
plot(svd, type = "scores", scores = 13:14)
plot(svd, type = "loadings", loadings = 1:10, coeff = 0.4)
plot(svd, type = "loadings", loadings = 11:20, coeff = 0.4)
plot(svd, type = "scores", scores = 13:14)
library(ggplot2)
plot(svd, type = "scores", scores = c(13, 15)) +
  aes(color = as.factor(G[, which.max(abs(svd$v[, 15]))])) +
  labs(color = "Geno")

ind.excl <- snp_indLRLDR(CHR, POS)
ind.keep <- snp_clumping(G, CHR, exclude = ind.excl, ncores = nb_cores())
svd2 <- big_randomSVD(G, snp_scaleBinom(), ind.col = ind.keep,
                      k = 20, ncores = nb_cores())
plot(svd2)
plot(svd2, type = "scores")
plot(svd2, type = "scores", scores = 3:4)
plot(svd2, type = "scores", scores = 13:14)
plot(svd2, type = "loadings", loadings = 1:10, coeff = 0.4)
plot(svd2, type = "loadings", loadings = 11:20, coeff = 0.4)

svd3 <- snp_autoSVD(G, CHR, POS, is.size.in.bp = TRUE,
                    k = 20, ncores = nb_cores())
plot(svd3, type = "loadings", loadings = 11:20, coeff = 0.4)
round(100 * cor(svd2$u, svd3$u)^2, 1)
