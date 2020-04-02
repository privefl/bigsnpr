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
# snp_subset(celiac, ind.row = which(log(dist) < 4), ind.col = which(CHR == 6),
#            backingfile = "../Dubois2010_data/celiac_chr6")

chr6 <- snp_attach("../Dubois2010_data/celiac_chr6.rds")
G <- chr6$genotypes$copy(code = c(0, 1, 2, 0, rep(NA, 252)))
dim(G)
big_counts(G, ind.col = 1:10)
CHR <- chr6$map$chromosome
POS <- chr6$map$physical.pos

POS2 <- snp_asGeneticPos(CHR, POS, dir = "tmp-data/")
plot(POS, POS2, pch = 20)

corr <- snp_cor(chr6$genotypes, infos.pos = POS2, size = 3 / 1000, ncores = 6)

# Simu phenotype
ind.HLA <- snp_indLRLDR(CHR, POS, LD.wiki34[12, ])
set.seed(1)
y2 <- snp_simuPheno(G, h2 = 0.1, M = 500)$pheno
y2 <- snp_simuPheno(G, h2 = 0.2, M = 1000, ind.possible = ind.HLA)$pheno


# GWAS
ind.gwas <- sample(nrow(G), 8e3)
gwas <- big_univLinReg(G, y2[ind.gwas], ind.train = ind.gwas)
beta_gwas <- gwas$estim
plot(gwas, type = "Manhattan")

ind.val <- setdiff(rows_along(G), ind.gwas)

# LDpred-inf
N <- length(ind.gwas)
df_beta <- data.frame(beta = gwas$estim, beta_se = gwas$std.err, n_eff = N)
beta_inf <- snp_ldpred2_inf(corr, df_beta)
pred <- big_prodVec(G, beta_inf, ind.row = ind.val)
plot(pred, y2[ind.val], pch = 20); abline(0, 1, col = "red", lwd = 3)
cor(pred, y2[ind.val])**2

# LDpred-grid
(ldsc <- snp_ldsc2(corr, df_beta, intercept = 1))
(h2_seq <- round(ldsc[["h2"]] + -2:2 * ldsc[["h2_se"]], 3))
(p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2))
(params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE)))

all_ldpred2 <- snp_ldpred2_grid(corr, df_beta, grid_param = params, ncores = 4)

pred2 <- big_prodMat(G, all_ldpred2, ind.row = ind.val)
params$r2 <- r2 <- drop(cor(pred2, y2[ind.val])**2)
round(100 * r2, 2)
round(100 * drop(cor(pred, y2[ind.val])**2), 2)

library(ggplot2)
ggplot(params, aes(x = p, y = r2, color = factor(h2))) +
  theme_bigstatsr() +
  geom_point() +
  geom_line() +
  scale_x_log10(breaks = 10^(-5:0), minor_breaks = params$p) +
  facet_wrap(~ sparse, labeller = label_both)

library(dplyr)
params %>%
  mutate(sparsity = colMeans(all_ldpred2 == 0)) %>%
  arrange(desc(r2)) %>%
  mutate_at(4:5, round, digits = 3) %>%
  slice(1:10)

# LDpred-auto
ldpred_auto <- snp_ldpred2_auto(corr, df_beta, verbose = TRUE)
pred3 <- big_prodVec(G, ldpred_auto[[1]]$beta_est, ind.row = ind.val)
round(100 * drop(cor(pred3, y2[ind.val])**2), 2)
