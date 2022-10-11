library(bigsnpr)

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
# Matrix::rankMatrix(corr, method = "qr")
# 18937  (out of 18941)

hist(log10(abs(corr@x)))

corr2 <- as_SFBM(corr)


set.seed(1)
# Simu phenotype
ind <- sample(ncol(G), ncol(G) * 0.2)
h2_1 <- 0.2
h2_2 <- 0.3
rho <- -0.4
rho_to_cov <- rho * sqrt(h2_1 * h2_2)
(cov <- matrix(c(h2_1, rho_to_cov, rho_to_cov, h2_2), 2))
beta <- MASS::mvrnorm(length(ind), rep(0, 2), cov / length(ind))

g <- big_prodMat(G, beta, ind.col = ind,
                 center = stats$center[ind],
                 scale = stats$scale[ind])
cov(g)

rhoE <- 0
rhoE_to_cov <- rhoE * sqrt((1 - h2_1) * (1 - h2_2))
covE <- matrix(c(1 - h2_1, rhoE_to_cov, rhoE_to_cov, 1 - h2_2), 2)
e <- MASS::mvrnorm(nrow(G), rep(0, 2), covE)
cov(y <- g + e)

# GWAS
ind.gwas <- sample(nrow(G), 7e3)
gwas <- big_univLinReg(G, y[ind.gwas, 1], ind.train = ind.gwas)
plot(gwas, type = "Manhattan")
hist(lpval <- -predict(gwas))

ind.val <- setdiff(rows_along(G), ind.gwas)

# LDSc reg
df_beta <- data.frame(beta = gwas$estim, beta_se = gwas$std.err,
                      n_eff = length(ind.gwas))
(ldsc <- snp_ldsc2(corr, df_beta))
ldsc_h2_est <- ldsc[["h2"]]


# LDpred2-auto
multi_auto <- snp_ldpred2_auto(corr2, df_beta, h2_init = ldsc_h2_est,
                               vec_p_init = seq_log(1e-4, 0.2, 8), ncores = 4,
                               burn_in = 300, num_iter = 200, report_step = 10)

(range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
(keep <- (range > (0.9 * quantile(range, 0.9))))

multi_auto_keep <- multi_auto[keep]

gamma <- rowMeans(sapply(multi_auto_keep, function(auto) auto$beta_est))
true_gamma <- rep(0, length(gamma)); true_gamma[ind] <- beta[, 1]
cor(gamma, true_gamma)

scale <- with(df_beta, sqrt(n_eff * beta_se^2 + beta^2))
all_bsamp <- lapply(multi_auto_keep, function(auto)
  as.matrix(sweep(auto$samp, 1, scale, '*')))

all_cov <- sapply(seq_along(all_bsamp), function(ic) {
  sapply(all_bsamp[-ic], function(other_bsamp)
    cov(all_bsamp[[ic]], other_bsamp))
})
hist(all_cov)
mean(all_cov)
abline(v = print(cov(gamma, true_gamma)), col = "red")

all_r2 <- sapply(seq_along(all_bsamp), function(ic) {
  one_mat <- do.call("cbind", all_bsamp[-ic])
  Rb <- apply(one_mat, 2, function(x) bigsparser::sp_prodVec(corr2, x))
  crossprod(all_bsamp[[ic]], Rb)
})
hist(all_r2)
pred <- big_prodVec(G, gamma, ind.row = ind.val)
abline(v = print(cor(pred, y[ind.val, 1]) ** 2), col = "red")
