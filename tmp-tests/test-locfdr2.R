df_beta <- readRDS("tmp-data/sumstats_tuto.rds")
corr    <- readRDS("tmp-data/corr_tuto.rds")
ld      <- readRDS("tmp-data/ld_scores_tuto.rds")

library(bigsnpr)
(p_seq  <- signif(seq_log(1, 1e-5, length.out = 10), 2))
(params <- expand.grid(h2 = 0.3, sparse = c(FALSE, TRUE), p = p_seq))

r <- with(df_beta, beta / sqrt(n_eff * beta_se^2 + beta^2))
pval <- with(df_beta, pchisq((beta / beta_se)^2, df = 1, lower.tail = FALSE))
z <- with(df_beta, beta / beta_se)

locfdr1 <- fdrtool::fdrtool(pval, statistic = "pvalue", plot = FALSE)$lfdr
locfdr2 <- fdrtool::fdrtool(r, statistic = "correlation", plot = FALSE)$lfdr
plot(locfdr1, locfdr2)
locfdr3 <- fdrtool::fdrtool(z, statistic = "normal", plot = FALSE)$lfdr
plot(locfdr1, locfdr3)
locfdr4 <- fdrtool::fdrtool(z, statistic = "normal",
                            cutoff.method = "locfdr", plot = FALSE)$lfdr
plot(locfdr1, locfdr4)



plot(pval, 1 - locfdr1, log = "x")
hist(pval)
hist(1 - locfdr1)
mean(1 - locfdr1)
