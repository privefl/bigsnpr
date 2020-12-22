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

corr <- runonce::save_run(
  snp_cor(chr6$genotypes, infos.pos = POS2, size = 3 / 1000, ncores = 6),
  file = "tmp-data/corr_chr6.rds"
)

# Simu phenotype
ind.HLA <- snp_indLRLDR(CHR, POS, LD.wiki34[12, ])
set.seed(1)
# simu <- snp_simuPheno(G, h2 = 0.5, M = 1000, ind.possible = ind.HLA)
simu <- snp_simuPheno(G, h2 = 0.2, M = 1000)
y2 <- simu$pheno


# GWAS
ind.gwas <- sample(nrow(G), 8e3)
gwas <- big_univLinReg(G, y2[ind.gwas], ind.train = ind.gwas)
plot(gwas, type = "Manhattan")
hist(lpval <- -predict(gwas))

ind.val <- setdiff(rows_along(G), ind.gwas)

# LDSc reg
df_beta <- data.frame(beta = gwas$estim, beta_se = gwas$std.err, n_eff = length(ind.gwas))
(ldsc <- snp_ldsc2(corr, df_beta))
h2_est <- ldsc[["h2"]]

THR <- 0.1
ind <- which(lpval > -log10(THR))
df_beta2 <- df_beta[ind, ]
z <- df_beta2$beta / df_beta2$beta_se

Z <- seq(0, 2 * max(abs(z)), length.out = 1e6)
thr <- sqrt(qchisq(THR, df = 1, lower.tail = FALSE))
Z2 <- Z + (dnorm(Z - thr) - dnorm(-Z - thr)) / (pnorm(Z - thr) + pnorm(-Z - thr))

df_beta3 <- df_beta2
knn <- nabor::knn(Z2, abs(z), k = 1)
new_z <- Z[drop(knn$nn.idx)] * sign(z)
plot(new_z, z); abline(0, 1, col = "red")
df_beta3$beta <- new_z * df_beta3$beta_se
plot(df_beta3$beta, df_beta2$beta); abline(0, 1, col = "red")


corr2 <- bigsparser::as_SFBM(as(corr[ind, ind], "dgCMatrix"))

(ldsc <- snp_ldsc2(corr[ind, ind], df_beta3))

# LDpred2-inf
beta_inf <- snp_ldpred2_inf(corr2, df_beta3, h2_est)
plot(beta_inf, df_beta$beta[ind])

# LDpred2-grid
(p_seq <- signif(seq_log(1e-4, 1, length.out = 17), 2))
(params <- expand.grid(p = p_seq, h2 = ldsc[["h2"]] * c(0.7, 1, 1.4), sparse = FALSE))
beta_grid <- snp_ldpred2_grid(corr2, df_beta3, params, ncores = 4)
pred_grid <- big_prodMat(G, beta_grid, ind.col = ind)
params$score <- big_univLinReg(as_FBM(pred_grid[ind.val, ]), y2[ind.val])$score
max(params$score)

library(ggplot2)
ggplot(params, aes(x = p, y = score, color = as.factor(h2))) +
  theme_bigstatsr() +
  geom_point() +
  geom_line() +
  scale_x_log10(breaks = 10^(-5:0), minor_breaks = params$p) +
  facet_wrap(~ sparse, labeller = label_both) +
  labs(y = "GLM Z-Score", color = "h2") +
  theme(legend.position = "top", panel.spacing = unit(1, "lines"))

# LDpred2-auto
auto <- snp_ldpred2_auto(corr2, df_beta3, h2_init = h2_est,
                         burn_in = 500, num_iter = 200, verbose = TRUE)
# with(auto[[1]], plot(beta_est, beta_est_sparse, pch = 20)); abline(0, 1, col = "red")
# mean(auto[[1]]$beta_est_sp == 0)  # 85%
# auto[[1]]$p_est



library(ggplot2)
postp <- auto[[1]]$postp_est
qplot(postp, geom = "density", fill = seq_along(postp) %in% simu$set, alpha = I(0.4)) +
  theme_bigstatsr() +
  scale_x_log10() +
  labs(fill = "Causal?")

1 /mean(1 / c(0.3, 0.3, 0.1))
1 /mean(1 / c(0.1, 0.1, 0.9))
log(mean(exp(c(0.3, 0.3, 0.1))))
log(mean(exp(c(0.1, 0.1, 0.9))))
