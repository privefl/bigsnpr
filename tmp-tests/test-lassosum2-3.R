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
stats <- big_scale()(G)

POS2 <- snp_asGeneticPos(CHR, POS, dir = "tmp-data/")
plot(POS, POS2, pch = 20)

corr <- runonce::save_run(
  snp_cor(chr6$genotypes, infos.pos = POS2, size = 3 / 1000, ncores = 6),
  file = "tmp-data/corr_chr6.rds"
)

# hist(log10(abs(corr@x)))

corr2 <- bigsparser::as_SFBM(as(corr, "dgCMatrix"))


set.seed(1)
# Simu phenotype
# ind <- sample(ncol(G), ncol(G) * 10^-runif(1, 0.5, 2.5))
ind <- sample(ncol(G), ncol(G) * 10^-runif(1, 0, 0.5))
h2_1 <- 0.2
h2_2 <- 0.2
rho <- -0.4
rho_to_cov <- rho * sqrt(h2_1 * h2_2)
(cov <- matrix(c(h2_1, rho_to_cov, rho_to_cov, h2_2), 2))
beta <- MASS::mvrnorm(length(ind), rep(0, 2), cov / length(ind))

g <- big_prodMat(G, beta, ind.col = ind,
                 center = stats$center[ind],
                 scale = stats$scale[ind])
cov(g)
e <- MASS::mvrnorm(nrow(G), rep(0, 2), matrix(c(1 - h2_1, 0, 0, 1 - h2_2), 2))
cov(y <- g + e)

# GWAS
ind.gwas <- sample(nrow(G), 8e3)
ind.val <- setdiff(rows_along(G), ind.gwas)
gwas <- big_univLinReg(G, y[ind.gwas, 1], ind.train = ind.gwas)
plot(gwas, type = "Manhattan")
hist(lpval <- -predict(gwas))

df_beta <- data.frame(beta = gwas$estim, beta_se = gwas$std.err,
                      n_eff = length(ind.gwas))

# Lassosum2
beta_lassosum2 <- snp_lassosum2(corr2, df_beta, ncores = 4)

params <- attr(beta_lassosum2, "grid_param")

pred_grid <- big_prodMat(G, beta_lassosum2, ind.row = ind.val)
params$score <- big_univLinReg(as_FBM(pred_grid), y[ind.val])$score
params

library(ggplot2)
ggplot(params, aes(x = lambda, y = score, color = as.factor(s))) +
  theme_bigstatsr() +
  geom_point() +
  geom_line() +
  scale_x_log10(breaks = 10^(-5:0)) +
  labs(y = "GLM Z-Score", color = "s") +
  theme(legend.position = "top", panel.spacing = unit(1, "lines"))

ggplot(params, aes(x = lambda, y = sparsity, color = as.factor(s))) +
  theme_bigstatsr() +
  geom_point() +
  geom_line() +
  scale_x_log10(breaks = 10^(-5:0)) +
  labs(y = "Sparsity", color = "s") +
  theme(legend.position = "top", panel.spacing = unit(1, "lines"))

ggplot(params, aes(x = lambda, y = num_iter, color = as.factor(s))) +
  theme_bigstatsr() +
  geom_point() +
  geom_line() +
  scale_x_log10(breaks = 10^(-5:0)) +
  labs(y = "Number of iterations", color = "s") +
  theme(legend.position = "top", panel.spacing = unit(1, "lines"))



# LDSc reg + LDpred2-auto
(ldsc <- snp_ldsc2(corr, df_beta))
ldsc_h2_est <- ldsc[["h2"]]
auto <- snp_ldpred2_auto(corr2, df_beta, h2_init = ldsc_h2_est,
                         burn_in = 300, num_iter = 200, report_step = 10,
                         verbose = TRUE)

pred <- big_prodVec(G, auto[[1]]$beta_est, ind.row = ind.val)
cor(pred, y[ind.val, 1]) ** 2
summary(lm(y[ind.val, 1] ~ pred))$coef["pred", "t value"]
max(params$score)
