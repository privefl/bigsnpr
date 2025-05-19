library(bigsnpr)
library(dplyr)

celiac <- snp_attach("../Dubois2010_data/FinnuncorrNLITUK1UK3hap300_QC_norel.rds")
dim(G <- celiac$genotypes)  # 15155 x 281122
CHR <- celiac$map$chromosome
POS <- celiac$map$physical.pos
y <- celiac$fam$affection - 1
sex <- celiac$fam$sex
NCORES <- nb_cores()
ind_chr6 <- which(CHR == 6)

pop <- snp_getSampleInfos(celiac, "../Dubois2010_data/FinnNLITUK1UK3.clusterv2")[[1]]
pop2 <- c("Netherlands", "Italy", "UK1", "UK2", "Finland")[pop]
table(pop2, exclude = NULL)
# Finland       Italy Netherlands         UK1         UK2
#    2436        1035        1623        3325        6736

big_counts(G, ind.col = 1:10)
G2 <- snp_fastImputeSimple(G, method = "mean2", ncores = NCORES)

obj.svd <- runonce::save_run(
  snp_autoSVD(G2, CHR, POS, ncores = NCORES),
  file = "tmp-data/celiac_svd.rds"
)
plot(obj.svd)
plot(obj.svd, type = "scores", scores = 1:6, coef = 0.6)
PC <- predict(obj.svd)[, 1:6]

pop_centers <- bigutilsr::geometric_median(PC, by_grp = pop2)

ldist_to_centers <- apply(pop_centers, 1, function(center) {
  log(rowSums(sweep(PC, 2, center, '-')^2))
})

ind_Fin <- which(ldist_to_centers[, "Finland"] < 7 & pop2 == "Finland")
ind_It <- which(ldist_to_centers[, "Italy"] < 6 & pop2 == "Italy")

corr_It <- runonce::save_run(file = "tmp-data/LD_It_chr6.rds", {
  POS2 <- snp_asGeneticPos(CHR[ind_chr6], POS[ind_chr6])
  corr0 <- snp_cor(G, ind.row = ind_It, ind.col = ind_chr6,
                   size = 3 / 1000, infos.pos = POS2, ncores = NCORES)
  as_SFBM(corr0, backingfile = "tmp-data/LD_It_chr6", compact = TRUE)
})

corr_Fin <- runonce::save_run(file = "tmp-data/LD_Fin_chr6.rds", {
  POS2 <- snp_asGeneticPos(CHR[ind_chr6], POS[ind_chr6])
  corr0 <- snp_cor(G, ind.row = ind_Fin, ind.col = ind_chr6,
                   size = 3 / 1000, infos.pos = POS2, ncores = NCORES)
  as_SFBM(corr0, backingfile = "tmp-data/LD_Fin_chr6", compact = TRUE)
})

ind_ItFin <- c(ind_It, ind_Fin)
corr_ItFin <- runonce::save_run(file = "tmp-data/LD_ItFin_chr6.rds", {
  POS2 <- snp_asGeneticPos(CHR[ind_chr6], POS[ind_chr6])
  corr0 <- snp_cor(G, ind.row = ind_ItFin, ind.col = ind_chr6,
                   size = 3 / 1000, infos.pos = POS2, ncores = NCORES)
  as_SFBM(corr0, backingfile = "tmp-data/LD_ItFin_chr6", compact = TRUE)
})

get_postp <- function(pop) {
  ind <- get(paste0("ind_", pop))
  gwas <- big_univLogReg(G2, y[ind], ind.train = ind, ind.col = ind_chr6,
                         covar.train = cbind(sex, PC)[ind, ])
  df_beta <- dplyr::transmute(gwas, beta = estim, beta_se = std.err, n_eff = length(ind))

  corr <- get(paste0("corr_", pop))
  ldpred_auto0 <- snp_ldpred2_auto(corr, df_beta, h2_init = 0.2, vec_p_init = 0.01,
                                   burn_in = 100, num_iter = 200, verbose = TRUE,
                                   allow_jump_sign = FALSE, shrink_corr = 0.95,
                                   use_MLE = FALSE)[[1]]
  ldpred_auto0$postp_est
}

postp_It <- get_postp("It")
summary(postp_It)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.003458 0.003690 0.004188 0.006615 0.005710 0.897038

postp_Fin <- get_postp("Fin")
summary(postp_Fin)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.006014 0.006613 0.007705 0.013554 0.010911 1.000000

postp_ItFin <- get_postp("ItFin")
summary(postp_ItFin)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.005426 0.005999 0.007040 0.014546 0.010238 1.000000


library(ggplot2)
library(dplyr)
tibble(postp_It, postp_Fin, postp_ItFin) %>%
  tidyr::pivot_longer(cols = 1:3) %>%
  ggplot(aes(value, fill = name)) +
  geom_density(alpha = 0.4) +
  scale_x_log10()


get_beta_hat <- function(ind) {
  gwas <- big_univLogReg(G2, y[ind], ind.train = ind, ind.col = ind_chr6,
                         covar.train = cbind(sex, PC)[ind, ])
  with(gwas, estim / sqrt(length(ind) * std.err^2 + estim^2))
}
beta_hat <- cbind(get_beta_hat(ind_It), get_beta_hat(ind_Fin))
plot(beta_hat)

var <- cbind(big_colstats(G2, ind.row = ind_It, ind.col = ind_chr6)$var,
             big_colstats(G2, ind.row = ind_Fin, ind.col = ind_chr6)$var)
plot(var)

mean_ld <- mean(
  bigsnpr:::ld_scores_sfbm(corr_It, ind_sub = cols_along(corr_It) - 1L, ncores = 2))

Rcpp::sourceCpp("tmp-tests/proto-ldpred2x.cpp")

ldpred_auto <- ldpred2x_gibbs(
  list_corr    = list(corr_It, corr_Fin),
  beta_hat     = beta_hat,
  n_vec        = matrix(c(length(ind_It), length(ind_Fin)),
                        nrow = ncol(corr_It), ncol = 2, byrow = TRUE),
  log_var      = log(var),
  ind_sub      = cols_along(corr_It) - 1L,
  p_init       = 0.01,
  h2_init      = 0.5,
  burn_in      = 100,
  num_iter     = 200,
  verbose      = TRUE,
  no_jump_sign = TRUE,
  shrink_corr  = 0.95,
  use_mle      = FALSE,
  alpha_bounds = c(-1.5, 0.5) + 1,
  mean_ld      = mean_ld
)

postp_x <- ldpred_auto$postp_est

plot(postp_x, postp_ItFin, pch = 20); abline(0, 1, col = "red", lwd = 2)
plot(tibble(postp_It, postp_Fin, postp_ItFin, postp_x))
plot(postp_x, postp_It, pch = 20); abline(0, 1, col = "red", lwd = 2)
plot(postp_x, postp_Fin, pch = 20); abline(0, 1, col = "red", lwd = 2)
