library(bigsnpr)

chr6 <- snp_attach("../Dubois2010_data/celiac_chr6.rds")
G <- chr6$genotypes$copy(code = c(0, 1, 2, rep(0, 253)))
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

corr2 <- as_SFBM(corr)


# Simu phenotype
set.seed(1)
simu <- snp_simuPheno(G, h2 = 0.05, M = ncol(G) * 0.005)

g <- big_prodVec(G, simu$allelic_effects, ind.col = simu$set)
var(g) # h2

y <- simu$pheno

# GWAS
ind.gwas <- sample(nrow(G), 7e3)
gwas <- big_univLinReg(G, y[ind.gwas], ind.train = ind.gwas)
plot(gwas, type = "Manhattan")
hist(lpval <- -predict(gwas))
N <- length(ind.gwas)

# Variance
df_beta <- data.frame(beta = gwas$estim, beta_se = gwas$std.err, n_eff = N)
scale <- with(df_beta, sqrt(n_eff * beta_se^2 + beta^2))
corr_hat <- df_beta$beta / scale
true_gamma <- rep(0, nrow(df_beta)); true_gamma[simu$set] <- simu$allelic_effects
corr_true <- true_gamma / scale
plot(bigsparser::sp_prodVec(corr2, corr_true), corr_hat); abline(0, 1, col = "red", lwd = 2)
eps <- bigsparser::sp_prodVec(corr2, corr_true) - corr_hat
cov(as.matrix(eps))
1 / N



# LDSc reg
(ldsc <- snp_ldsc2(corr, df_beta))


# LDpred2-auto
set.seed(1)
Niter <- 600
multi_auto <- snp_ldpred2_auto(corr2, df_beta, h2_init = ldsc[["h2"]],
                               vec_p_init = seq_log(1e-4, 0.2, 20), ncores = 4,
                               allow_jump_sign = FALSE, shrink_corr = 0.95,
                               burn_in = 500, num_iter = Niter, report_step = 2)

signif(range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))), 2)
(keep <- (range > (0.95 * quantile(range, 0.95))))


## Rhat (convergence diagnosis)

# h2
(all_h2_est <- sapply(multi_auto, function(auto) auto$h2_est))
hist(all_h2 <- sapply(multi_auto, function(auto) tail(auto$path_h2_est, Niter)))
rstan::Rhat(all_h2)
rstan::ess_bulk(all_h2)

(Rhat_h2 <- sapply(multi_auto, function(auto)
  rstan::Rhat(tail(auto$path_h2_est, Niter))))
plot(Rhat_h2, all_h2_est); abline(v = 1.1, col = "red")
rstan::ess_bulk(all_h2[, which(Rhat_h2 < 1.1)])

cbind(Rhat_h2, all_h2_est)
plot(multi_auto[[3]]$path_h2_est, main = paste("Rhat =", round(Rhat_h2[3], 3)))
split <- rstan:::z_scale(rstan:::split_chains(tail(multi_auto[[2]]$path_h2_est, Niter)))
hist(split[, 1]); hist(split[, 2], add = TRUE, col = scales::alpha("blue", 0.2))


(all_p_est <- sapply(multi_auto, function(auto) auto$p_est))
hist(all_p <- sapply(multi_auto, function(auto) tail(auto$path_p_est, Niter)))
rstan::Rhat(all_p)
rstan::ess_bulk(all_p)

(Rhat_p <- sapply(multi_auto, function(auto)
  rstan::Rhat(tail(auto$path_p_est, Niter))))
plot(Rhat_p, all_p_est); abline(v = 1.1, col = "red")
rstan::ess_bulk(all_p[, which(Rhat_p < 1.1)])
plot(multi_auto[[1]]$path_p_est, main = paste("Rhat =", round(Rhat_p[1], 3)))
plot(multi_auto[[8]]$path_p_est, main = paste("Rhat =", round(Rhat_p[8], 3)))
plot(multi_auto[[20]]$path_p_est, main = paste("Rhat =", round(Rhat_p[20], 3)))

plot(Rhat_h2, Rhat_p)

# (Rhat_alpha <- sapply(multi_auto, function(auto)
#   rstan::Rhat(tail(auto$path_alpha_est, Niter))))
# all_alpha <- sapply(multi_auto, function(auto) tail(auto$path_alpha_est, Niter))
# hist(all_alpha[, which(Rhat_alpha < 1.1)])
# (all_alpha_est <- sapply(multi_auto, function(auto) auto$alpha_est))
# rstan::ess_bulk(all_alpha[, which(Rhat_alpha < 1.1)])



(stat2 <- sapply(multi_auto, function(auto) {

  x <- tail(auto$path_h2_est, Niter)
  split <- rstan:::z_scale(rstan:::split_chains(x))
  # hist(split[, 1]); hist(split[, 2], add = TRUE, col = scales::alpha("blue", 0.2))
  # rstan::Rhat(x)

  stat1 <- goftest::cvm.test(split[, 1], "pnorm")$statistic
  stat2 <- goftest::cvm.test(split[, 2], "pnorm")$statistic
  # these seem to always have the same value; let's take the max anyway
  max(stat1, stat2)
}))

plot(Rhat_h2, stat2)

(stat3 <- sapply(multi_auto, function(auto) {
  x <- tail(auto$path_p_est, Niter)
  split <- rstan:::z_scale(rstan:::split_chains(x))
  # hist(split[, 1]); hist(split[, 2], add = TRUE, col = scales::alpha("blue", 0.2))
  rstan::Rhat(x)

  stat1 <- goftest::cvm.test(split[, 1], "pnorm")$statistic
  stat2 <- goftest::cvm.test(split[, 2], "pnorm")$statistic
  # these seem to always have the same value; let's take the max anyway
  max(stat1, stat2)
}))

plot(Rhat_p, stat3)


(stat4 <- sapply(multi_auto, function(auto) {
  x <- tail(auto$path_h2_est, Niter)
  r <- rank(x, na.last = "keep", ties.method = "average")
  N <- floor(length(x) / 2)
  M <- length(x) - N
  # Equation (6) of 10.1214/aoms/1177704477
  t <- (N / M * crossprod(sort(r[1:N]) - 1:N * (N + M) / N) +
          M / N * crossprod(sort(r[1:M + N]) - 1:M * (N + M) / M)) / (N + M)^2
  t[[1]]
}))
plot(stat4, stat2)

(stat5 <- sapply(multi_auto, function(auto) {
  x <- tail(auto$path_p_est, Niter)
  r <- rank(x, na.last = "keep", ties.method = "average")
  N <- M <- Niter / 2
  t <- (N / M * crossprod(sort(r[1:N]) - 1:N * (N + M) / N) +
          M / N * crossprod(sort(r[1:M + N]) - 1:M * (N + M) / M)) / (N + M)^2
  t[[1]]
}))
plot(stat5, stat3)

## Test for the null distribution

size <- 500
nsim <- 5000
X <- matrix(rnorm(size * nsim), size)

null_pval0 <- apply(X, 2, function(x) {
  goftest::cvm.test(x, "pnorm")$p.value
})
hist(null_pval0)

size <- 5000
N <- floor(size / 2)
M <- size - N
ind1 <- 1:N; ind2 <- 1:M + N
exp_rank1 <- 1:N * (N + M) / N
exp_rank2 <- 1:M * (N + M) / M

null_stat0 <- replicate(n = 500e3, {
  r <- sample(size)
  (N / M * crossprod(sort(r[ind1]) - exp_rank1) +
      M / N * crossprod(sort(r[ind2]) - exp_rank2)) / (N + M)^2
})
hist(null_stat0)
mean(null_stat0 - (1 + 1 / size) / 6)
# hist(log(null_stat0))
exp(bigutilsr::tukey_mc_up(log(null_stat0), coef = 1.5))
# for 5000: 1.357353
# For 500: 1.348421 / 1.351163 / 1.359333
# For 50: 1.383114
exp(bigutilsr::tukey_mc_up(log(null_stat0), coef = 2))
# For 5000: 2.53935
# For 500: 2.522322 / 2.526508 / 2.544412
# For 50: 2.585694
quantile(null_stat0, 0.999)
# For 5000: 1.167146
# For 500: 1.171352
# For 50: 1.126802


null_pval <- apply(X, 2, function(x) {
  split <- rstan:::z_scale(rstan:::split_chains(x))
  # hist(split[, 1]); hist(split[, 2], add = TRUE, col = scales::alpha("blue", 0.2))
  # rstan::Rhat(x)
  pval1 <- goftest::cvm.test(split[, 1], "pnorm")$p.value
  pval2 <- goftest::cvm.test(split[, 2], "pnorm")$p.value
  # these seem to always have the same value; let's take the max anyway
  min(pval1, pval2)
})
hist(null_pval)

rint <- function(x) {
  # x[] <- qnorm(tail(dplyr::percent_rank(c(-Inf, Inf, x)), -2))
  r <- rank(x, na.last = "keep", ties.method = "average")
  p <- (r - 3/8) / (sum(!is.na(r)) + 1/4)
  x[] <- qnorm(p)
  x
}
rint(print(matrix(c(1:3, NA), 2)))
rint(print(matrix(c(1, 1, 1, NA), 2)))

null_pval2 <- apply(X, 2, function(x) {
  split <- rint(rstan:::split_chains(x))
  # hist(split[, 1]); hist(split[, 2], add = TRUE, col = scales::alpha("blue", 0.2))
  rstan::Rhat(x)
  pval1 <- goftest::cvm.test(split[, 1], "pnorm")$p.value
  pval2 <- goftest::cvm.test(split[, 2], "pnorm")$p.value
  # these seem to always have the same value; let's take the max anyway
  min(pval1, pval2)
})
hist(null_pval2)


# can detect a trend
x <- 1:300
goftest::cvm.test(rint(x), "pnorm")$stat
split <- matrix(rint(x), ncol = 2)
goftest::cvm.test(split[, 1], "pnorm")$stat
goftest::cvm.test(split[, 2], "pnorm")$stat
rstan::Rhat(x)


# almost constant
x <- rep(0, 1000); x[c(50, 572)] <- 1
plot(x)
rstan::Rhat(x)
split <- rstan:::z_scale(rstan:::split_chains(x))
goftest::cvm.test(split[, 1], "pnorm")$statistic
goftest::cvm.test(split[, 2], "pnorm")$statistic
r <- rank(x, na.last = "keep", ties.method = "average")
N <- M <- length(x) / 2
t <- (N / M * crossprod(sort(r[1:N]) - 1:N * (N + M) / N) +
        M / N * crossprod(sort(r[1:M + N]) - 1:M * (N + M) / M)) / (N + M)^2
t[[1]]



## Partial likelihood
chi2 <- with(df_beta, (beta / beta_se)^2)
ind_keep <- snp_clumping(G, infos.chr = CHR, infos.pos = POS, thr.r2 = 0.2, S = chi2)
ind <- ind_keep[order(chi2[ind_keep], decreasing = TRUE)[1:20]]
chi2[ind]
corr2[ind, ind]

sapply(multi_auto, function(auto) {
  corr_imp <- auto$corr_est
  diff <- (corr_hat - corr_imp)[ind]
  lik <- -diff %*% Matrix::solve(corr2[ind, ind] / N, diff)
  lik[1]
})


all_lik <- sapply(multi_auto, function(auto) {
  beta <- auto$sample_beta[ind, ]
  R <- corr2[ind, ind]
  corr_imp <- R %*% beta
  corr_hat_sub <- corr_hat[ind]
  apply(corr_imp, 2, function(x) {
    diff <- corr_hat_sub - x
    lik <- -diff %*% Matrix::solve(R, diff)
    lik[1]
  })
})

mystat <- function(x) {
  split <- rstan:::z_scale(rstan:::split_chains(x))
  stat1 <- goftest::cvm.test(split[, 1], "pnorm")$statistic
  stat2 <- goftest::cvm.test(split[, 2], "pnorm")$statistic
  # these seem to always have the same value; let's take the max anyway
  max(stat1, stat2)
}
plot(apply(all_lik, 2, rstan::Rhat), apply(all_lik, 2, mystat))

hist(all_lik)
hist(exp(all_lik))
hist(all_lik[, 20], 50, main = paste("Rhat =", round(rstan::Rhat(all_lik[, 20]), 3)))
hist(all_lik[, 20], 50, main = paste("My Stat =", round(mystat(all_lik[, 20]), 3)))
rstan::Rhat(all_lik[, 20])
rstan::Rhat(all_lik[, 19])

plot(multi_auto[[20]]$path_p_est)

ind.val <- setdiff(rows_along(G), ind.gwas)
all_r2_test <- sapply(multi_auto, function(auto) {
  pred <- big_prodVec(G, auto$beta_est, ind.row = ind.val)
  cor(pred, y[ind.val])^2
})
plot(colMeans(all_lik), all_r2_test)
