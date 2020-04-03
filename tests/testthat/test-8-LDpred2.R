################################################################################

context("LDPRED2")

################################################################################

test_that("LDpred2 works", {

  skip_if(is_cran)
  skip_if_offline("dropbox.com")

  load(url("https://www.dropbox.com/s/c13uygnjh6yh7vf/to-test-ldpred2.RData?raw=1"))

  # LDpred2-inf
  beta_inf <- snp_ldpred2_inf(corr, df_beta)
  r_inf <- cor(beta_inf, true_beta)
  expect_gt(r_inf, 0.4)

  # Naive PRS
  expect_lt(cor(df_beta$beta, true_beta), r_inf)

  # LDpred2-gibbs
  p_seq <- signif(seq_log(1e-3, 1, length.out = 7), 1)
  ldsc <- snp_ldsc2(corr, df_beta, intercept = 1, blocks = NULL)
  params <- expand.grid(p = p_seq, h2 = ldsc[["h2"]], sparse = c(FALSE, TRUE))
  expect_equal(dim(params), c(14, 3))
  beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = 2)
  expect_gt(max(cor(beta_grid, true_beta)), r_inf)
  expect_gt(min(cor(beta_grid, true_beta)), 0.3)

  # LDpred2-auto
  beta_auto <- snp_ldpred2_auto(corr, df_beta, burn_in = 500, num_iter = 200)
  expect_gt(cor(beta_auto$beta_est, true_beta), r_inf)

  # Errors
  expect_error(snp_ldpred2_inf(corr, df_beta[-1]),
               "'df_beta' should have element 'beta'.")
  expect_error(snp_ldpred2_grid(corr, df_beta[-1], params),
               "'df_beta' should have element 'beta'.")
  expect_error(snp_ldpred2_auto(corr, df_beta[-1]),
               "'df_beta' should have element 'beta'.")

  # OpenMP
  all_beta <- replicate(
    n = 10, simplify = FALSE,
    snp_ldpred2_grid(corr, df_beta, params[c(6, 13), ], ncores = 2))
  lapply(all_beta[-1], function(.) expect_equal(., all_beta[[1]]))
})

################################################################################
