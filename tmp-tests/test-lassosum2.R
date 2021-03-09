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
ind <- sample(ncol(G), ncol(G) * 10^-runif(1, 0.5, 2.5))
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

# LDSc reg
(ldsc <- snp_ldsc2(corr, df_beta))
ldsc_h2_est <- ldsc[["h2"]]


# LDpred2-auto
auto <- snp_ldpred2_auto(corr2, df_beta, h2_init = ldsc_h2_est,
                         burn_in = 300, num_iter = 200, report_step = 10,
                         verbose = TRUE)

pred <- big_prodVec(G, auto[[1]]$beta_est, ind.row = ind.val)
cor(pred, y[ind.val, 1]) ** 2



# Lassosum2
N <- df_beta$n_eff
scale <- sqrt(N * df_beta$beta_se^2 + df_beta$beta^2)
beta_hat <- df_beta$beta / scale

Rcpp::sourceCpp('tmp-tests/test-lassosum2.cpp')
beta_lassosum2 <- lassosum2(
  corr = corr2,
  beta_hat = beta_hat,
  beta_init = 0 * beta_hat,
  order = seq_along(beta_hat) - 1L,
  lambda = 0.01,
  s = 0.9
)
mean(beta_lassosum2 == 0)

pred_lassosum2 <- big_prodVec(G, beta_lassosum2 * scale, ind.row = ind.val)
cor(pred_lassosum2, y[ind.val, 1]) ** 2

plot(beta_lassosum2 * scale, auto[[1]]$beta_est)


beta_lassosum0 <- lassosum2(
  corr = corr2,
  beta_hat = beta_hat,
  beta_init = 0 * beta_hat,
  order = seq_along(beta_hat) - 1L,
  lambda = 0.01,
  s = 1
)

pred_lassosum0 <- big_prodVec(G, beta_lassosum0 * scale, ind.row = ind.val)
cor(pred_lassosum0, y[ind.val, 1]) ** 2

plot(beta_lassosum0 * scale, auto[[1]]$beta_est)





N <- df_beta$n_eff
scale <- sqrt(N * df_beta$beta_se^2 + df_beta$beta^2)
beta_hat <- df_beta$beta / scale

lambda0 <- max(abs(beta_hat))
params <- expand.grid(
  lambda = seq_log(0.01 * lambda0, lambda0, 10),
  s = c(0.5, 0.8, 0.9, 0.95, 1)
)

library(furrr)
plan("multisession", workers = 5)

all_res <- future_pmap_dfr(params, function(lambda, s) {

  Rcpp::sourceCpp('tmp-tests/test-lassosum2.cpp')

  time <- system.time(
    beta_lassosum <- lassosum2(
      corr = corr2,
      beta_hat = beta_hat,
      beta_init = 0 * beta_hat,
      order = seq_along(beta_hat) - 1L,
      lambda = lambda,
      s = s
    )
  )[3]

  pred_lassosum <- big_prodVec(G, beta_lassosum * scale, ind.row = ind.val)
  r2 <- cor(pred_lassosum, y[ind.val, 1]) ** 2

  tibble::tibble(lambda, s, time, sparsity = mean(beta_lassosum == 0), r2)
}, .progress = TRUE)
# saveRDS(all_res, "tmp-data/all_res_to_plot.rds")


# time is larger for smaller n and smaller lambda
# all_res <- readRDS("tmp-data/all_res_to_plot.rds")
all_res
max(all_res$r2)

library(ggplot2)
ggplot(all_res, aes(x = lambda, y = r2, color = as.factor(s))) +
  bigstatsr::theme_bigstatsr() +
  geom_point() +
  geom_line() +
  scale_x_log10(breaks = 10^(-5:0)) +
  labs(y = "r2", color = "s") +
  theme(legend.position = "top")

library(ggplot2)
ggplot(all_res, aes(x = lambda, y = time, color = as.factor(s))) +
  bigstatsr::theme_bigstatsr() +
  geom_point() +
  geom_line() +
  scale_x_log10(breaks = 10^(-5:0)) +
  labs(y = "Time", color = "s") +
  theme(legend.position = "top")


plot(beta_lassosum3 * scale, auto[[1]]$beta_est)

# warm stats does not seem to work
# need to choose sequence of lambdas and s
# and see if do not have to always compute crossprods
# use different lambdas
