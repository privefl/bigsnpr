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
y2 <- snp_simuPheno(G, h2 = 0.2, M = 5000)$pheno


# GWAS
ind.gwas <- sample(nrow(G), 8e3)
gwas <- big_univLinReg(G, y2[ind.gwas], ind.train = ind.gwas)
plot(gwas, type = "Manhattan")

ind.val <- setdiff(rows_along(G), ind.gwas)

# LDSc reg
df_beta <- data.frame(beta = gwas$estim, beta_se = gwas$std.err, n_eff = length(ind.gwas))
(ldsc <- snp_ldsc2(corr, df_beta))
h2_est <- ldsc[["h2"]]


# LDpred-inf
beta_inf <- snp_ldpred2_inf(corr, df_beta)
pred_inf <- big_prodVec(G, beta_inf, ind.row = ind.val)
mad(pred_inf)
plot(pred_inf, y2[ind.val], pch = 20); abline(0, 1, col = "red", lwd = 3)
cor(pred_inf, y2[ind.val])**2

# LDpred-grid
(p_seq <- signif(seq_log(1e-4, 1, length.out = 9), 2))
(params <- expand.grid(p = p_seq, h2 = h2_est, sparse = c(FALSE, TRUE)))
beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = 6)
pred_grid <- big_prodMat(G, beta_grid, ind.row = ind.val)
params$r <- apply(pred_grid, 2, cor, y = y[ind.val])
params$coef <- apply(pred_grid, 2, function(x) {
  glm(y[ind.val] ~ x)$coef[2]
})
params$mad <- apply(pred_grid, 2, mad)
library(ggplot2)
plot_grid(
  ggplot(params, aes(x = p, y = coef, color = as.factor(h2))) +
    theme_bigstatsr() +
    geom_point() +
    geom_line() +
    scale_x_log10(breaks = 10^(-5:0), minor_breaks = params$p) +
    facet_wrap(~ sparse, labeller = label_both) +
    labs(y = "coef", color = "h2") +
    theme(legend.position = "top", panel.spacing = unit(1, "lines")),

  ggplot(params, aes(x = p, y = r, color = as.factor(h2))) +
    theme_bigstatsr() +
    geom_point() +
    geom_line() +
    scale_x_log10(breaks = 10^(-5:0), minor_breaks = params$p) +
    facet_wrap(~ sparse, labeller = label_both) +
    labs(y = "corr", color = "h2") +
    theme(legend.position = "top", panel.spacing = unit(1, "lines")),

  ggplot(params, aes(x = p, y = mad, color = as.factor(h2))) +
    theme_bigstatsr() +
    geom_point() +
    geom_line() +
    scale_x_log10(breaks = 10^(-5:0), minor_breaks = params$p) +
    facet_wrap(~ sparse, labeller = label_both) +
    labs(y = "MAD", color = "h2") +
    theme(legend.position = "top", panel.spacing = unit(1, "lines")),

  nrow = 2
)

# LDpred-auto
auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                         burn_in = 200, num_iter = 200, verbose = TRUE)
pred_auto <- big_prodVec(G, auto$beta_est, ind.row = ind.val)
mad(pred_auto)
plot(pred_auto, y2[ind.val], pch = 20); abline(0, 1, col = "red", lwd = 3)
cor(pred_auto, y2[ind.val])**2

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

logit <- function(p) log(p / (1 - p))
curve(log(x))

path_log_p <- -log(auto$path_p_est[201:400])
median(path_log_p) / mad(path_log_p)
path_log_h2 <- -log(auto$path_h2_est[201:400])
median(path_log_h2) / mad(path_log_h2)

library(future.apply)
plan(multisession(workers = 4))
multi_auto <- future_lapply(c(0.001, 0.01, 0.1, 0.9), function(p) {
  auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                           p_init = p, burn_in = 1000, num_iter = 1000)
})

sapply(multi_auto, function(auto) {
  beta_auto <- auto$beta_est
  pred_auto <- big_prodVec(G, beta_auto, ind.row = ind.val)
  mad(pred_auto)
})

library(ggplot2)
auto <- multi_auto[[1]]
for (auto in multi_auto) {
  print(plot_grid(
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
  ))
}

