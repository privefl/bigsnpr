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
ind <- sample(ncol(G), ncol(G) * 10^-runif(1, 0.5, 2.5))
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
auto <- snp_ldpred2_auto(corr2, df_beta, h2_init = ldsc_h2_est,
                         burn_in = 300, num_iter = 200, report_step = 10,
                         verbose = FALSE)

scale <- with(df_beta, sqrt(n_eff * beta_se^2 + beta^2))
pred_samp <- big_prodVec(G, auto[[1]]$sample_beta[, 1] * scale, ind.row = ind.val)
var(pred_samp) # ~h2_1
plot(g[ind.val, 1], pred_samp); abline(0, 1, col = "red")
lm(pred_samp ~ g[ind.val, 1])$coef

pred_samp2 <- big_prodVec(G, auto[[1]]$sample_beta[, 2] * scale, ind.row = ind.val)
var(pred_samp2) # ~h2_1
lm(pred_samp2 ~ g[ind.val, 1])$coef

pred_samp3 <- big_prodVec(G, rowMeans(auto[[1]]$sample_beta) * scale, ind.row = ind.val)
var(pred_samp3)
all_alpha <- apply(auto[[1]]$sample_beta, 2, function(x) {
  pred_x <- big_prodVec(G, x * scale, ind.row = ind.val)
  lm(pred_x ~ g[ind.val, 1])$coef[[2]]
})
mean(all_alpha * h2_1)

mod1 <- lm(pred_samp ~ g[ind.val, 1]); alpha1 <- mod1$coef[[2]]
mod2 <- lm(pred_samp2 ~ g[ind.val, 1]); alpha2 <- mod2$coef[[2]]
c(h2_1 * alpha1 * alpha2 + cov(mod1$resid, mod2$resid),
  cov(pred_samp, pred_samp2))
# c(cov(mod1$resid, mod2$resid), (1 - alpha1) * (1 - alpha2) * h2_1^2)

c(cov(pred_samp,  y[ind.val, 1]), h2_1 * alpha1)
c(cov(pred_samp2, y[ind.val, 1]), h2_1 * alpha2)
c(cor(pred_samp, y[ind.val, 1])^2, h2_1 * alpha1^2)
c(cor(pred_samp2, y[ind.val, 1])^2, h2_1 * alpha2^2)

pred <- big_prodVec(G, auto[[1]]$beta_est, ind.row = ind.val)
c(crossprod(auto[[1]]$sample_beta[, 2],
            bigsparser::sp_prodVec(corr2, auto[[1]]$sample_beta[, 1])),
  cor(pred, y[ind.val, 1]) ** 2)
