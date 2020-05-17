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

corr <- snp_cor(chr6$genotypes, infos.pos = POS2, size = 3 / 1000,
                ncores = 6, alpha = 0.9)
object.size(corr) / 1024**2  # 0.3 -> 51 Mb / 0.9 -> 79 / 1 -> 84
str(corr)
median(Matrix::colSums(corr != 0))

ld <- Matrix::colSums(corr ** 2)
plot(ld, pch = 20)
plot(POS, bigutilsr::rollmean(ld, 100))
hist(S <- bigutilsr::rollmean(ld, 100), "FD")
abline(v = (q <- bigutilsr::tukey_mc_up(S)), col = "red")

# Simu phenotype
ind.HLA <- snp_indLRLDR(CHR, POS, LD.wiki34[12, ])
set.seed(1)
y2 <- snp_simuPheno(G, h2 = 0.5, M = 1000, ind.possible = ind.HLA)$pheno
# y2 <- snp_simuPheno(G, h2 = 0.02, M = 5000)$pheno


# GWAS
ind.gwas <- sample(nrow(G), 8e3)
gwas <- big_univLinReg(G, y2[ind.gwas], ind.train = ind.gwas)
plot(gwas, type = "Manhattan")

ind.val <- setdiff(rows_along(G), ind.gwas)

# LDSc reg
df_beta <- data.frame(beta = gwas$estim, beta_se = gwas$std.err, n_eff = length(ind.gwas))
(ldsc <- snp_ldsc2(corr, df_beta))
h2_est <- ldsc[["h2"]]

m <- ncol(corr)
N <- df_beta$n_eff
sd <- df_beta$beta_se * sqrt(N)
beta_hat <- df_beta$beta / sd
system.time(
  beta_inf1 <- as.vector(Matrix::solve(
    corr + Matrix::Diagonal(m, m / (h2_est * N)), beta_hat))
)
cor(big_prodVec(G, beta_inf1 * sd, ind.row = ind.val), y2[ind.val])**2

corr2 <- bigsparser::as_SFBM(as(corr, "dgCMatrix"))
system.time(
  beta_inf2 <- bigsparser::sp_solve_sym(corr2, beta_hat,
                                        add_to_diag = m / (h2_est * N))
)
cor(big_prodVec(G, beta_inf2 * sd, ind.row = ind.val), y2[ind.val])**2
plot(beta_inf1, beta_inf2)
crossprod(beta_inf2, bigsparser::sp_prodVec(corr2, beta_inf2))

# LDpred-inf
beta_inf <- snp_ldpred2_inf(corr2, df_beta, h2 = h2_est)
pred_inf <- big_prodVec(G, beta_inf, ind.row = ind.val)
mad(pred_inf)
plot(pred_inf, y2[ind.val], pch = 20); abline(0, 1, col = "red", lwd = 3)
cor(pred_inf, y2[ind.val])**2

# LDpred-grid
(p_seq <- signif(seq_log(1e-4, 1, length.out = 9), 2))
(params <- expand.grid(p = p_seq, h2 = h2_est, sparse = c(FALSE, TRUE)))
beta_grid <- snp_ldpred2_grid(corr2, df_beta, params, ncores = 6)
pred_grid <- big_prodMat(G, beta_grid, ind.row = ind.val)
params$r <- apply(pred_grid, 2, cor, y = y2[ind.val])

params[c("coef", "score")] <-
  big_univLinReg(as_FBM(pred_grid), y2[ind.val])[c("estim", "score")]

library(ggplot2)
ggplot(params, aes(x = p, y = score, color = as.factor(h2))) +
  theme_bigstatsr() +
  geom_point() +
  geom_line() +
  scale_x_log10(breaks = 10^(-5:0), minor_breaks = params$p) +
  facet_wrap(~ sparse, labeller = label_both) +
  labs(y = "GLM Z-Score", color = "h2") +
  theme(legend.position = "top", panel.spacing = unit(1, "lines"))

params$h2_est2 <- apply(beta_grid, 2, function(beta) {
  crossprod(beta / sd, bigsparser::sp_prodVec(corr2, beta / sd))
})

ggplot(params, aes(x = p, y = h2_est, color = as.factor(h2))) +
  theme_bigstatsr() +
  geom_point() +
  geom_line() +
  scale_x_log10(breaks = 10^(-5:0), minor_breaks = params$p) +
  facet_wrap(~ sparse, labeller = label_both) +
  labs(y = "h2_est", color = "h2") +
  theme(legend.position = "top", panel.spacing = unit(1, "lines"))

ggplot(params, aes(x = p, y = r, color = as.factor(h2))) +
  theme_bigstatsr() +
  geom_point() +
  geom_line() +
  scale_x_log10(breaks = 10^(-5:0), minor_breaks = params$p) +
  facet_wrap(~ sparse, labeller = label_both) +
  labs(y = "Correlation", color = "h2") +
  theme(legend.position = "top", panel.spacing = unit(1, "lines"))


# LDpred-auto
auto <- snp_ldpred2_auto(corr2, df_beta, h2_init = h2_est,
                         burn_in = 200, num_iter = 200, verbose = TRUE)
pred_auto <- big_prodVec(G, auto$beta_est, ind.row = ind.val)
var(pred_auto)
plot(pred_auto, y2[ind.val], pch = 20); abline(0, 1, col = "red", lwd = 3)
cor(pred_auto, y2[ind.val])**2

crossprod(auto$beta_est / sd, bigsparser::sp_prodVec(corr2, auto$beta_est / sd))

library(ggplot2)
plot_grid(
  qplot(y = auto$path_p_est) +
    theme_bigstatsr() +
    geom_hline(yintercept = auto$p_est, col = "blue") +
    scale_y_log10() +
    labs(y = "p"),
  qplot(y = auto$path_h2_est) +
    theme_bigstatsr() +
    geom_hline(yintercept = auto$h2_est, col = "blue") +
    labs(y = "h2"),
  ncol = 1, align = "hv"
)

library(future.apply)
plan(multisession(workers = 6))

library(dplyr)
res <- tibble(p_init = seq_log(1e-3, 0.9, 12))

multi_auto <- future_lapply(res$p_init, function(p) {
  auto <- snp_ldpred2_auto(corr2, df_beta, h2_init = h2_est,
                           p_init = p, burn_in = 500, num_iter = 200)
})

res$cp <- sapply(multi_auto, function(auto) {
  beta_auto <- auto$beta_est / sd
  crossprod(beta_auto, bigsparser::sp_prodVec(corr2, beta_auto))
}) %>% print()

library(ggplot2)
auto <- multi_auto[[1]]
for (auto in multi_auto) {
  print(auto[2:3])
  print(plot_grid(
    qplot(y = auto$path_p_est) +
      theme_bigstatsr() +
      geom_hline(yintercept = auto$p_est, col = "blue") +
      scale_y_log10() +
      labs(y = "p"),
    qplot(y = auto$path_h2_est) +
      theme_bigstatsr() +
      geom_hline(yintercept = auto$h2_est, col = "blue") +
      scale_y_log10() +
      labs(y = "h2"),
    ncol = 1, align = "hv"
  ))
}

betas_auto <- sapply(multi_auto, function(auto) auto$beta_est)
cor(betas_auto)

d <- dist(t(betas_auto))
hc <- hclust(d)
plot(hc)

res$dist_to_others <- colMeans(as.matrix(d**2))

preds_auto <- big_prodMat(G, betas_auto, ind.row = ind.val)
(res$perf <- apply(preds_auto, 2, cor, y = y2[ind.val])**2)

res$h2_est <- purrr::map_dbl(multi_auto, "h2_est")
res$p_est <- purrr::map_dbl(multi_auto, "p_est")

plot(res[-1])

summary(lm(perf ~ ., data = res[-1]))

ggplot(res, aes(cp, perf, color = dist_to_others)) +
  geom_point() +
  theme_bigstatsr() +
  geom_abline() +
  labs(color = NULL) +
  scale_colour_viridis_c(direction = 1)
plot(colMeans(as.matrix(d**2)), cor_auto)


ggplot(res, aes(p_init, p_est, color = dist_to_others)) +
  geom_point() +
  theme_bigstatsr() +
  geom_abline() +
  labs(color = NULL) +
  scale_colour_viridis_c(direction = 1)


h2_auto <- sapply(multi_auto, function(auto)
  auto$beta_est * bigsparser::sp_prodVec(corr2, auto$beta_est))
d <- dist(t(h2_auto))
hc <- hclust(d)
plot(hc)
colMeans(as.matrix(d**2))

plot(h2_auto[, 1], h2_auto[, 2])

pca <- snp_autoSVD(G, CHR, POS, ncores = 6)
ind.keep <- attr(pca, "subset")
obj.svd <- big_SVD(as_FBM(t(h2_auto)), big_scale(), ind.col = ind.keep)
# obj.svd <- big_SVD(as_FBM(t(h2_auto)), ind.col = ind.keep)
plot(obj.svd)
plot(obj.svd, type = "scores")
plot(obj.svd, type = "loadings")

dist_pc <- sqrt(rowSums(obj.svd$u[, 1:2]**2))
plot(dist_pc, cor_auto)

sapply(multi_auto, function(auto)
  stats::acf(tail(auto$path_p_est, 200), lag.max = 1, plot = FALSE)$acf)
sapply(multi_auto, function(auto)
  stats::acf(tail(auto$path_p_est, 200), lag.max = 10, plot = FALSE)$acf)
res$acf <- colMeans(sapply(multi_auto, function(auto)
  stats::acf(tail(auto$path_p_est, 200), lag.max = 10, plot = FALSE)$acf))
sapply(multi_auto, function(auto)
  stats::acf(tail(auto$path_h2_est, 200), lag.max = 1, plot = FALSE)$acf)
res$acf2 <- sapply(multi_auto, function(auto)
  stats::acf(tail(auto$path_p_est, 200), lag.max = 10, plot = FALSE)$acf)[11, ]

ggplot(res, aes(acf, perf, color = dist_to_others)) +
  geom_point() +
  theme_bigstatsr() +
  geom_abline() +
  labs(color = NULL) +
  scale_colour_viridis_c(direction = 1)
