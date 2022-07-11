library(bigsnpr)
options(max.print = 500)

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

POS2 <- snp_asGeneticPos(CHR, POS, dir = "tmp-data")
plot(POS, POS2, pch = 20)

corr <- runonce::save_run(
  snp_cor(chr6$genotypes, infos.pos = POS2, size = 3 / 1000, ncores = 6),
  file = "tmp-data/corr_chr6.rds"
)

corr2 <- bigsparser::as_SFBM(corr, compact = TRUE)

ind <- 5000 + -2:1
round(corr[ind, ind], 3)

# Simu phenotype
# set.seed(1)
simu <- snp_simuPheno(G, h2 = 0.05, M = length(ind), ind.possible = ind)
simu$effects
y2 <- simu$pheno


# GWAS
ind.gwas <- sample(nrow(G), 8e3)
gwas <- big_univLinReg(G, y2[ind.gwas], ind.train = ind.gwas)
# plot(gwas, type = "Manhattan")

ind.val <- setdiff(rows_along(G), ind.gwas)

# LDSc reg
df_beta <- data.frame(beta = gwas$estim, beta_se = gwas$std.err, n_eff = length(ind.gwas))
(ldsc <- snp_ldsc2(corr, df_beta))
h2_est <- ldsc[["h2"]]


# LDpred-auto
auto <- snp_ldpred2_auto(corr2, df_beta, h2_init = h2_est,
                          burn_in = 100, num_iter = 100, report_step = 5,
                          allow_jump_sign = FALSE, shrink_corr = 0.95,
                          verbose = TRUE)[[1]]
# Overall: 0.000142106 // 0.0509896

pred2 <- big_prodVec(G, auto$beta_est, ind.row = ind.val)
cor(pred2, y2[ind.val])^2

plot(auto$postp_est)

scale <- with(df_beta, sqrt(n_eff * beta_se^2 + beta^2))
beta_hat <- df_beta$beta / scale
plot(beta_hat, auto$corr_est); abline(0, 1, col = "red", lwd = 2)

diff <- beta_hat - auto$corr_est
hist(diff)
newZ <- diff * scale / df_beta$beta_se
hist(newZ, probability = TRUE); curve(dnorm(x), col = "red", add = TRUE)
bsamp <- auto$sample_beta
cbind(true = simu$effects, bsamp[ind, ])



beta_hat <- with(df_beta[ind, ], beta /  sqrt(n_eff * beta_se^2 + beta^2))

simu$effects
bigsparser::sp_solve_sym(
  as_SFBM(corr[ind, ind]), beta_hat, add_to_diag = 0.005)


ind2 <- 4900:5100
beta_hat2 <- with(df_beta[ind2, ], beta /  sqrt(n_eff * beta_se^2 + beta^2))

betas_inf <- bigsparser::sp_solve_sym(
  as_SFBM(corr[ind2, ind2]), beta_hat2, add_to_diag = 0.005)
betas_inf[match(ind, ind2)]
simu$effects

beta_inf <- snp_ldpred2_inf(corr2, df_beta, h2 = 0.05 * 10)
beta_inf[ind] * with(df_beta[ind, ], sqrt(n_eff * beta_se^2 + beta^2))
simu$effects
