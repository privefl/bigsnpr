library(bigsnpr)

# celiac <- snp_attach("../Dubois2010_data/FinnuncorrNLITUK1UK3hap300_QC_norel.rds")
# G <- celiac$genotypes$copy(code = c(0, 1, 2, 0, rep(NA, 252)))
# CHR <- celiac$map$chromosome
# POS <- celiac$map$physical.pos
#
# obj.svd <- snp_autoSVD(G, CHR, POS, ncores = 4, thr.r2 = 0.1)
# plot(obj.svd)
#
# dist <- bigutilsr::dist_ogk(obj.svd$u)
# hist(log(dist))
#
# snp_subset(celiac, ind.row = which(log(dist) < 4), ind.col = which(CHR == 22),
#            backingfile = "../Dubois2010_data/celiac_chr22")

chr22 <- snp_attach("../Dubois2010_data/celiac_chr22.rds")
G <- chr22$genotypes$copy(code = c(0, 1, 2, 0, rep(NA, 252)))
dim(G)  # 11402  4945
CHR <- chr22$map$chromosome
POS <- chr22$map$physical.pos
POS2 <- snp_asGeneticPos(CHR, POS, dir = "tmp-data")

corr <- snp_cor(chr22$genotypes, infos.pos = POS2, size = 3 / 1000, ncores = 4)
object.size(corr) / 1024^2  # 17 Mb

# Simu phenotype
set.seed(1)
simu <- snp_simuPheno(G, 0.3, 1000)
sd <- big_scale()(G, rows_along(G), simu$set)$scale
true_beta <- rep(0, ncol(G))
true_beta[simu$set] <- simu$effects / sd

# GWAS
y <- simu$pheno
gwas <- big_univLinReg(G, y)
df_beta <- data.frame(
  beta    = gwas$estim,
  beta_se = gwas$std.err,
  n_eff   = nrow(G))

y2 <- (y > qnorm(0.8)) + 0L
gwas2 <- big_univLogReg(G, y2)
df_beta2 <- data.frame(
  beta    = gwas2$estim,
  beta_se = gwas2$std.err,
  n_eff   = 4 / (1 / sum(y2 == 0) + 1 / sum(y2 == 1)))

save(corr, df_beta, df_beta2, true_beta,
     file = "tmp-data/to-test-ldpred2.RData")
