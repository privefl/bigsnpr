df_beta <- readRDS("tmp-data/sumstats_tuto.rds")
corr    <- readRDS("tmp-data/corr_tuto.rds")
ld      <- readRDS("tmp-data/ld_scores_tuto.rds")

pval <- with(df_beta, pchisq((beta / beta_se)^2, df = 1, lower.tail = FALSE))

hist(pval)

cat("Step 1... determine cutoff point\n")

x0 = fdrtool::fndr.cutoff(pval, "pvalue")

cat("Step 2... estimate parameters of null distribution and eta0\n")

eta0 <- fdrtool::censored.fit(x = pval, cutoff = x0, statistic = "pvalue")[1, 3]

cat("Step 3... compute p-values and estimate empirical PDF/CDF\n")

ee <- fdrtool:::ecdf.pval(pval, eta0 = eta0)
ee(seq(0, 1, by = 0.02))
g.pval <- fdrtool::grenander(ee)
plot(g.pval)

plot(g.pval$x.knots, g.pval$f.knots, log = "xy")
f.pval = approxfun(g.pval$x.knots, g.pval$f.knots, method = "constant",
                   rule = 2)

plot(f.pval, log = "y", ylim = c(eta0, 2))


f.pval2 = function(x) {
  exp(spline(log(g.pval$x.knots[-1]), log(g.pval$f.knots[-1]), xout = log(x),
             method = "hyman")$y)
}
plot(f.pval2, log = "y")

fdr.pval = function(p) {
  p[p == .Machine$double.eps] = 0
  pmin(eta0/f.pval(p), 1)
}
fdr.pval2 = function(p) {
  p[p == .Machine$double.eps] = 0
  pmin(eta0/f.pval2(p), 1)
}

cat("Step 4... compute local fdr\n")

lfdr <- fdr.pval(pval)
lfdr2 <- fdr.pval2(pval)

plot(lfdr, lfdr2, log = "xy")

plot(pval, 1 - lfdr, log = "x")
points(pval, 1 - lfdr2, pch = 20, col = "red")

ord <- order(pval)
print(1 - lfdr[ord], digits = 20)

mean(lfdr2 == 1)
mean(1 - lfdr2)  # 0.01142471


auto <- bigsnpr::snp_ldpred2_auto(corr, df_beta, h2_init = 0.3, verbose = TRUE)[[1]]
auto$p_est # 0.01210376

library(ggplot2)
qplot(auto$postp_est, 1 - lfdr2, color = ld) +
  theme_bw() +
  # scale_x_log10() +
  scale_color_viridis_c(trans = "log", direction = -1)

roll <- function(x) bigutilsr::rollmean(x, size = 10)

qplot(roll(auto$postp_est), roll(1 - lfdr2), color = ld) +
  theme_bw() +
  # scale_x_log10() +
  scale_color_viridis_c(trans = "log", direction = -1)


qplot(roll(auto$postp_est), roll(1 - lfdr2) / sqrt(ld), color = ld) +
  theme_bw() +
  # scale_x_log10() +
  scale_color_viridis_c(trans = "log", direction = -1) +
  geom_abline(color = "red", linetype = 2)


beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = 0.3)

library(bigsnpr)
# Read from bed/bim/fam, it generates .bk and .rds files.
snp_readBed("tmp-data/public-data3.bed")
obj.bigSNP <- snp_attach("tmp-data/public-data3.rds")
G   <- obj.bigSNP$genotypes
y   <- obj.bigSNP$fam$affection

pred_inf <- big_prodVec(G, beta_inf, ind.col = df_beta[["_NUM_ID_"]])
cor(pred_inf, y)  # 0.2638186


p <- pmax(1 - lfdr2, 1e-5)

N <- df_beta$n_eff
scale <- sqrt(N * df_beta$beta_se^2 + df_beta$beta^2)
beta_hat <- df_beta$beta / scale

beta_inf2 <- scale * bigsparser::sp_solve_sym(
  corr, beta_hat, add_to_diag = sum(p) / (0.3 * N * p))

pred_inf2 <- big_prodVec(G, beta_inf2, ind.col = df_beta[["_NUM_ID_"]])
cor(pred_inf2, y)  # 0.39 with min_p = 1e-5  //  0.44 with min_p = 0.1


p <- (1 - lfdr2) / ld
h2 <- 0.3
m <- ncol(corr)
coeff <- N * h2 / m
beta_blup <- beta_inf / scale
C1 <- (coeff + 1) / (coeff + p)
C2 <- coeff / (coeff + 1)
V_nc <- C2 / N * (1 - C2 * h2)
V_c <- C2 * h2 / (m * p) + V_nc
d_beta_ncaus <- dnorm(beta_blup, sd = sqrt(V_nc))
d_beta_caus <- dnorm(beta_blup, sd = sqrt(V_c))
d_caus_beta <- d_beta_caus * p / (d_beta_caus * p + d_beta_ncaus * (1 - p))
beta_jms <- beta_inf * C1 * d_caus_beta
plot(beta_inf, beta_jms)

pred_inf3 <- big_prodVec(G, beta_jms, ind.col = df_beta[["_NUM_ID_"]])
cor(pred_inf3, y)  # 0.417 with p -- 0.432 with p/ld



multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = 0.3,
                               vec_p_init = seq_log(1e-4, 0.5, length.out = 15),
                               allow_jump_sign = FALSE, shrink_corr = 0.95,
                               ncores = nb_cores())

beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)
pred_auto <- big_prodMat(G, beta_auto, ind.col = df_beta[["_NUM_ID_"]])
sc <- apply(pred_auto, 2, sd)
keep <- abs(sc - median(sc)) < 3 * mad(sc)
final_pred_auto <- rowMeans(pred_auto[, keep])
cor(final_pred_auto, y)  # 0.4863071
