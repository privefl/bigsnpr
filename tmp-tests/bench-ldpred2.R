df_beta <- readRDS("tmp-data/sumstats_tuto.rds")
corr    <- readRDS("tmp-data/corr_tuto.rds")

library(bigsnpr)
(p_seq  <- signif(seq_log(1, 1e-5, length.out = 10), 2))
(params <- expand.grid(h2 = 0.3, sparse = c(FALSE, TRUE), p = p_seq))

system.time(
  beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = 4)
)
# 27-30 / 7.8-8.5

# cor(beta_grid)


SEQ_P <- signif(seq_log(1e-5, 1, length.out = 6), 2)
all_time <- sapply(SEQ_P, function(p) {
  print(p)
  params <- expand.grid(p = rep(p, 4), h2 = 0.3, sparse = c(FALSE, TRUE))
  system.time(
    beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = 4)
  )
})
plot(SEQ_P, print(all_time[3, ]), log = "x")
# Without UnrolOpt: 3.92 3.89 4.00 3.97 4.24 6.61
# With UnrolOpt: 3.25 3.34 3.19 3.34 3.51 5.71
# With test for 0: 3.27 3.25 3.32 3.31 3.67 5.88
