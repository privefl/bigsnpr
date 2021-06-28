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

corr2 <- bigsparser::as_SFBM(corr)


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


#### lassosum2 ####
beta_lassosum2 <- snp_lassosum2(corr2, df_beta,
                                ncores = 4,
                                lambda.min.ratio = 0.05,
                                nlambda = 10,
                                maxiter = 100)

params <- attr(beta_lassosum2, "grid_param")

pred_grid <- big_prodMat(G, beta_lassosum2, ind.row = ind.val)
params$score <- apply(pred_grid, 2, cor, y = y[ind.val])
dplyr::arrange(params, dplyr::desc(score))

library(ggplot2)
ggplot(params, aes(x = lambda, y = score, color = as.factor(delta))) +
  theme_bigstatsr() +
  geom_point() +
  geom_line() +
  scale_x_log10(breaks = 10^(-5:0)) +
  labs(y = "Correlation", color = "delta") +
  theme(legend.position = "top", panel.spacing = unit(1, "lines"))


#### PRS-CS ####

N <- df_beta$n_eff
scale <- sqrt(N * df_beta$beta_se^2 + df_beta$beta^2)

Rcpp::sourceCpp('tmp-tests/test-prscs.cpp')

system.time(
  test <- prscs_gibbs_one(
    sfbm = corr2$address,
    beta_hat = df_beta$beta / scale,
    n_vec = N,
    a = 1,
    b = 0.5,
    phi = 0.01,
    max_psi = 1,
    burn_in = 100,
    num_iter = 100
  )
) # 4 sec -> 5 sec (when sampling from gig)

pred_prscs <- big_prodVec(G, test * scale, ind.row = ind.val)
cor(pred_prscs, y[ind.val])

cov(g)

params2 <- expand.grid(a = c(0.5, 1, 1.5),
                       b = c(0, 0.5, 1),
                       phi = seq_log(1e-5, 1, 6))

test2 <- prscs_gibbs(
  corr2,
  beta_hat = df_beta$beta / scale,
  # order = seq_along(scale) - 1L,
  n_vec = N,
  a = params2$a,
  b = params2$b,
  phi = params2$phi,
  max_psi = 1,
  burn_in = 100,
  num_iter = 100,
  ncores = 1  # returns NAs when using parallelism, WTF?
)

pred_prscs2 <- big_prodMat(G, sweep(test2, 1, scale, '*'), ind.row = ind.val)
params2$score <- cor(pred_prscs2, y[ind.val])
dplyr::arrange(params2, dplyr::desc(score))

library(ggplot2)
ggplot(params2, aes(x = phi, y = score, color = paste(a, b, sep = " & "))) +
  theme_bigstatsr() +
  geom_point() +
  geom_line() +
  scale_x_log10(breaks = 10^(-5:0)) +
  labs(y = "Correlation", color = "a & b") +
  theme(legend.position = "top", panel.spacing = unit(1, "lines"))
