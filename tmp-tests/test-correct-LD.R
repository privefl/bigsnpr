library(bigsnpr)
options(max.print = 500)

correct <- function(corr, n, pval_thr) {

  r <- corr@x

  z <- sqrt(r^2 * (n - 2) / (1 - r^2))

  max_z <- 0.9 * sqrt((n - 2) / (1 - 0.9^2))

  thr_Z <- sqrt(stats::qchisq(pval_thr, df = 1, lower.tail = FALSE))
  Z <- seq(0, max_z, length.out = 1e5)
  Z2 <- Z + (stats::dnorm(Z - thr_Z) - stats::dnorm(-Z - thr_Z)) /
    (stats::pnorm(Z - thr_Z) + stats::pnorm(-Z - thr_Z))

  # plot(Z2, Z, log = "xy"); abline(0, 1, col = "red", lwd = 2)

  knn <- bigutilsr::knn_parallel(Z2, as.matrix(z), k = 1, ncores = 1)
  new_z <- Z[drop(knn$nn.idx)]
  corr@x <- ifelse(
    abs(r) > 0.9, r, ifelse(
      z < thr_Z, 0, sign(r) * new_z / sqrt(new_z^2 + n - 2)))

  Matrix::drop0(corr)
}

chr6 <- snp_attach("../Dubois2010_data/celiac_chr6.rds")
G <- chr6$genotypes$copy(code = c(0, 1, 2, 0, rep(NA, 252)))
dim(G)
big_counts(G, ind.col = 1:10)
CHR <- chr6$map$chromosome
POS <- chr6$map$physical.pos

corr <- runonce::save_run({
  POS2 <- snp_asGeneticPos(CHR, POS, dir = "tmp-data")
  plot(POS, POS2, pch = 20)
  snp_cor(chr6$genotypes, infos.pos = POS2, size = 3 / 1000, ncores = 6)
}, file = "tmp-data/corr_chr6.rds")


# Simu phenotype
set.seed(1)
simu <- snp_simuPheno(G, h2 = 0.1, M = 100)
simu$effects
y2 <- simu$pheno


# GWAS
ind.gwas <- sample(nrow(G), 8e3)
gwas <- big_univLinReg(G, y2[ind.gwas], ind.train = ind.gwas)
# plot(gwas, type = "Manhattan")

ind.val <- setdiff(rows_along(G), ind.gwas)

# LDSc reg
df_beta <- data.frame(beta = gwas$estim, beta_se = gwas$std.err, n_eff = length(ind.gwas))
(ldsc <- snp_ldsc2(corr, df_beta))


SEQ <- 10^-(0:10)
ALL_RES <- sapply(SEQ, function(pval_thr) {
  print(pval_thr)
  corr_thr <- correct(corr, n, pval_thr)
  round(cbind(corr[1:5, 1:5], corr_thr[1:5, 1:5]), 3)
  corr2 <- bigsparser::as_SFBM(corr_thr, compact = TRUE)
  time <- system.time(
    auto <- snp_ldpred2_auto(corr2, df_beta, h2_init = ldsc[["h2"]])[[1]]
  )[[3]]
  pred <- big_prodVec(G, auto$beta_est, ind.row = ind.val)
  c(time, cor(pred, y2[ind.val]))
})
plot(SEQ, ALL_RES[1, ], log = "x")
plot(SEQ, ALL_RES[2, ], log = "x")
