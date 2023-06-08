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
# set.seed(1)
simu <- snp_simuPheno(G, h2 = 0.03, M = ncol(G) * 0.01)

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
                               burn_in = 500, num_iter = Niter, report_step = 10)

signif(range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))), 2)
(keep <- (range > (0.95 * quantile(range, 0.95))))


# Rhat (convergence diagnosis)
all_h2 <- sapply(multi_auto, function(auto) tail(auto$path_h2_est, Niter))
rstan::Rhat(all_h2)
rstan::ess_bulk(all_h2)
hist(all_h2)

all_p <- sapply(multi_auto, function(auto) tail(auto$path_p_est, Niter))
rstan::Rhat(all_p)
rstan::ess_bulk(all_p)
hist(all_p)

(all_h2_est <- sapply(multi_auto, function(auto) auto$h2_est))
(all_p_est <- sapply(multi_auto, function(auto) auto$p_est))

(Rhat_h2 <- sapply(multi_auto, function(auto)
  rstan::Rhat(tail(auto$path_h2_est, Niter))))
rstan::ess_bulk(all_h2[, which(Rhat_h2 < 1.1)])

cbind(Rhat_h2, all_h2_est)
plot(multi_auto[[3]]$path_h2_est, main = paste("Rhat =", round(Rhat_h2[3], 3)))
split <- rstan:::z_scale(rstan:::split_chains(tail(multi_auto[[3]]$path_h2_est, Niter)))
hist(split)

(Rhat_p <- sapply(multi_auto, function(auto)
  rstan::Rhat(tail(auto$path_p_est, Niter))))
rstan::ess_bulk(all_p[, which(Rhat_p < 1.1)])

(Rhat_alpha <- sapply(multi_auto, function(auto)
  rstan::Rhat(tail(auto$path_alpha_est, Niter))))
all_alpha <- sapply(multi_auto, function(auto) tail(auto$path_alpha_est, Niter))
hist(all_alpha[, which(Rhat_alpha < 1.1)])
(all_alpha_est <- sapply(multi_auto, function(auto) auto$alpha_est))
rstan::ess_bulk(all_alpha[, which(Rhat_alpha < 1.1)])

plot(Rhat_h2, all_h2_est); abline(v = 1.1, col = "red")
plot(Rhat_p, all_p_est); abline(v = 1.1, col = "red")


kept_h2 <- rstan:::split_chains(all_h2)

all_h2_z <- rstan:::z_scale(kept_h2)
stats <- apply(all_h2_z, 2, function(x)
  `if`(anyNA(x), NA, goftest::cvm.test(x, "pnorm", mean = 0, sd = 1)$statistic))
print(stats)
ind <- which.max(stats)
kept_h2[, ind] <- NA
# [1]  98.237645  54.591009   6.616759  51.166368 120.221137   4.998669  52.238196  28.414644
# [9]  54.246587   1.151838 141.923347  14.108193  30.741331  26.884741   9.569849  14.344620
# [17]  39.817221   5.798867  64.362072  23.629123  80.173934  25.418182  19.950714 101.950336
# [25] 150.921418   9.318483   4.623866  44.887964  57.936784  23.513504 120.591376 108.322916
# [33]  30.984441  35.029280   8.045592  26.652380   6.503355   8.319606   3.322661  17.663176
# [1] 107.569173  47.634339   6.923594  57.154283 127.830069   5.717949  58.646128  23.707176
# [9]  47.032140   2.018210 150.478405  10.499751  25.855588  21.696751   6.881799  11.329155
# [17]  45.763306   5.631244  71.798080  20.881396  88.234146  21.453067  16.550714 111.758669
# [25]         NA   7.882803   4.442386  39.121637  50.973472  18.990776 131.170645  98.375613
# [33]  26.150286  30.081815   7.219605  30.774292   4.510370   7.011654   2.047209  18.369961

hist(rstan:::split_chains(all_h2)[, 5])

