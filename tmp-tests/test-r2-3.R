library(bigsnpr)

chr6 <- snp_attach("../Dubois2010_data/celiac_chr6.rds")
G <- chr6$genotypes$copy(code = c(0, 1, 2, rep(0, 253)))
dim(G)
big_counts(G, ind.col = 1:10)
CHR <- chr6$map$chromosome
POS <- chr6$map$physical.pos
stats <- big_scale()(G)

POS2 <- snp_asGeneticPos(CHR, POS, dir = "tmp-data")
plot(POS, POS2, pch = 20)

corr <- runonce::save_run(
  snp_cor(chr6$genotypes, infos.pos = POS2, size = 3 / 1000, ncores = 6),
  file = "tmp-data/corr_chr6.rds"
)

# Matrix::rankMatrix(corr, method = "qr")
# 18937  (out of 18941)

corr2 <- as_SFBM(corr)


# Simu phenotype
# set.seed(1)
simu <- snp_simuPheno(G, h2 = 0.01, M = ncol(G) * 0.005)

g <- big_prodVec(G, simu$allelic_effects, ind.col = simu$set)
var(g) # h2

y <- simu$pheno

# GWAS
ind.gwas <- sample(nrow(G), 7e3); ind.val <- setdiff(rows_along(G), ind.gwas)
gwas <- big_univLinReg(G, y[ind.gwas], ind.train = ind.gwas)
plot(gwas, type = "Manhattan")
hist(lpval <- -predict(gwas))
N <- length(ind.gwas)
df_beta <- data.frame(beta = gwas$estim, beta_se = gwas$std.err, n_eff = N)

scale <- with(df_beta, sqrt(n_eff * beta_se^2 + beta^2))
corr_hat <- df_beta$beta / scale


# LDSc reg
(ldsc <- snp_ldsc2(corr, df_beta))
ldsc_h2_est <- ldsc[["h2"]]


# LDpred2-auto
multi_auto <- snp_ldpred2_auto(corr2, df_beta, h2_init = ldsc_h2_est,
                               vec_p_init = seq_log(1e-4, 0.2, 20), ncores = 4,
                               burn_in = 300, num_iter = 200, report_step = 10)

(range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
(keep <- (range > (0.95 * quantile(range, 0.95))))

multi_auto_keep <- multi_auto[keep]

gamma <- rowMeans(sapply(multi_auto_keep, function(auto) auto$beta_est))
cor(gamma, true_gamma)
# Variance
true_gamma <- rep(0, length(gamma)); true_gamma[simu$set] <- simu$allelic_effects
corr_true <- true_gamma / scale
plot(bigsparser::sp_prodVec(corr2, corr_true), corr_hat); abline(0, 1, col = "red", lwd = 2)
eps <- bigsparser::sp_prodVec(corr2, corr_true) - corr_hat
cov(as.matrix(eps))
1 / N


all_bsamp <- lapply(multi_auto_keep, function(auto) as.matrix(auto$samp))

all_gsamp <- lapply(multi_auto_keep, function(auto)
  as.matrix(sweep(auto$samp, 1, scale, '*')))

all_cov <- sapply(seq_along(all_gsamp), function(ic) {
  sapply(all_gsamp[-ic], function(other_gsamp)
    cov(all_gsamp[[ic]], other_gsamp))
})
hist(all_cov)
mean(all_cov)
abline(v = print(cov(gamma, true_gamma)), col = "red")

all_r2 <- sapply(seq_along(all_bsamp), function(ic) {
  one_mat <- do.call("cbind", all_bsamp[-ic])
  Rb <- apply(one_mat, 2, function(x) bigsparser::sp_prodVec(corr2, x))
  crossprod(all_bsamp[[ic]], Rb)
})
hist(all_r2, "FD")
pred <- big_prodVec(G, gamma, ind.row = ind.val)
abline(v = print(cor(pred, y[ind.val]) ** 2), col = "red")

beta <- gamma / scale
bRb <- crossprod(beta, corr %*% beta)[1]
hist(all_r2^2 / bRb, "FD")
abline(v = print(cor(pred, y[ind.val]) ** 2), col = "red")

c(var(pred), bRb)
# I could use bRb as a criterion to select chains maybe? No, cf. below

all_bRb <- sapply(multi_auto, function(auto) {
  beta <- auto$beta_est / scale
  crossprod(beta, corr %*% beta)[1]
})
plot(all_bRb)
plot(all_bRb, range)  # no real corr


all_r2_test <- sapply(multi_auto, function(auto) {
  pred <- big_prodVec(G, auto$beta_est, ind.row = ind.val)
  cor(pred, y[ind.val])^2
})
plot(all_bRb, all_r2_test) # no real corr
plot(range, all_r2_test) # good positive corr


all_r2_2 <- sapply(multi_auto, function(auto) {
  bsamp <- auto$sample_beta
  Rb <- apply(bsamp, 2, function(x) bigsparser::sp_prodVec(corr2, x))
  b1Rb2 <- crossprod(bsamp, Rb)
  beta <- auto$beta_est / scale
  bRb <- crossprod(beta, corr %*% beta)[1]
  mean(as(Matrix::tril(b1Rb2, k = -2), "sparseMatrix")@x^2 / bRb)
})
plot(all_r2_2, all_r2_test)


cor(cbind(all_bRb, range, all_r2_2, all_r2_test), method = "spearman")

chi2 <- with(df_beta, (beta / beta_se)^2)
ind_keep <- snp_clumping(G, infos.chr = CHR, infos.pos = POS, thr.r2 = 0.2, S = chi2)
ind <- ind_keep[order(chi2[ind_keep], decreasing = TRUE)[1:1000]]
chi2[ind]
corr2[ind, ind]

all_lik <- sapply(multi_auto, function(auto) {
  corr_imp <- auto$corr_est
  diff <- (corr_hat - corr_imp)[ind]
  lik <- -diff %*% Matrix::solve(corr2[ind, ind] / N, diff)
  lik[1]
})
plot(range, all_r2_test)
plot(all_lik, all_r2_test)
cor(cbind(all_bRb, range, all_r2_2, all_lik, all_r2_test), method = "spearman")

all_lik2 <- sapply(multi_auto, function(auto) {
  corr_imp <- auto$corr_est
  sig <- corr2[ind, ind] / N
  k <- Matrix::rankMatrix(sig)
  det <- Matrix::determinant(sig, logarithm = TRUE)$modulus
  diff <- (corr_hat - corr_imp)[ind]
  lik <- -0.5 * (det + diff %*% Matrix::solve(sig, diff) + k * log(2 * pi))
  exp(lik[1])
})

plot(all_lik2, all_r2_test)

cor(cbind(all_bRb, range, all_r2_2, all_lik, all_lik2, all_r2_test), method = "spearman")

(h2_seq <- round(ldsc_h2_est * c(0.3, 0.7, 1, 1.4), 4))
(p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2))
(params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE)))
beta_grid <- snp_ldpred2_grid(corr2, df_beta, grid_param = params, ncores = 4)

all_lik_grid <- apply(beta_grid, 2, function(beta) {
  corr_imp <- beta / scale
  diff <- (corr_hat - corr_imp)[ind]
  lik <- -diff %*% Matrix::solve(corr2[ind, ind] / N, diff)
  lik[1]
})

all_r2_test_grid <- apply(beta_grid, 2, function(beta) {
  pred <- big_prodVec(G, beta, ind.row = ind.val)
  cor(pred, y[ind.val])^2
})

plot(all_lik_grid, all_r2_test_grid, col = params$sparse + 1)
