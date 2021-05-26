################################################################################

context("LDPRED2")

################################################################################

test_that("sp_colSumsSq_sym() works", {

  replicate(100, {

    N <- 300
    spmat <- Matrix::rsparsematrix(N, N, 0.1, symmetric = TRUE)
    expect_equal(sp_colSumsSq_sym(spmat@p, spmat@i, spmat@x),
                 Matrix::colSums(spmat^2))
  })

})

################################################################################

test_that("LDpred2 works", {

  skip_if(is_cran)
  skip_if_offline("dropbox.com")

  load(url("https://www.dropbox.com/s/c13uygnjh6yh7vf/to-test-ldpred2.RData?raw=1"))

  # LD score regression
  ldsc <- snp_ldsc2(corr, df_beta)
  corr <- as_SFBM(corr)

  # LDpred2-inf
  beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = ldsc[["h2"]])
  r_inf <- cor(beta_inf, true_beta)
  expect_gt(r_inf, 0.4)

  # Naive PRS
  expect_lt(cor(df_beta$beta, true_beta), r_inf)

  # LDpred2-gibbs
  p_seq <- signif(seq_log(1e-3, 1, length.out = 7), 1)
  params <- expand.grid(p = p_seq, h2 = ldsc[["h2"]] / 5, sparse = c(FALSE, TRUE))
  expect_equal(dim(params), c(14, 3))
  beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = 2)
  expect_gt(max(cor(beta_grid, true_beta), na.rm = TRUE), 0.4)

  # LDpred2-uncertainty
  expect_error(snp_ldpred2_grid(corr, df_beta, params, return_sampling_betas = TRUE),
               "Only one set of parameters is allowed")
  beta_sample <- snp_ldpred2_grid(corr, df_beta, params[2, ], num_iter = 200,
                                  return_sampling_betas = TRUE)
  expect_equal(dim(beta_sample), c(nrow(df_beta), 200))
  if (sd(rowMeans(beta_sample)) < 1) {
    ind <- which.max(cor(rowMeans(beta_sample), beta_grid[, 1:7]))
    expect_true(ind %in% 1:3)
    # not exactly the same, but should be close:
    expect_gt(cor(rowMeans(beta_sample), beta_grid[, ind]), 0.7)
  }

  # LDpred2-auto
  beta_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = ldsc[["h2"]],
                                burn_in = 200, num_iter = 200, sparse = TRUE)
  expect_length(beta_auto, 1)
  mod <- beta_auto[[1]]
  expect_null(dim(mod$beta_est))
  expect_gt(cor(mod$beta_est, true_beta), 0.4)
  expect_gt(cor(mod$beta_est_sparse, true_beta), 0.4)
  expect_lt(mean(mod$beta_est == 0), 0.001)
  expect_gt(mean(mod$beta_est_sparse == 0), 0.1)
  expect_equal(mod$h2_init, ldsc[["h2"]])
  expect_equal(mod$p_init, 0.1)
  expect_equal(mod$h2_est, mean(tail(mod$path_h2_est, 200)))
  expect_equal(mod$p_est,  mean(tail(mod$path_p_est,  200)))
  expect_length(mod$postp_est, length(mod$beta_est))
  expect_null(dim(mod$postp_est))
  expect_equal(mean(mod$postp_est), mod$p_est, tolerance = 0.01)
  expect_equal(dim(mod$sample_beta), c(ncol(corr), 0))

  # Sampling betas
  beta_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = ldsc[["h2"]],
                                burn_in = 200, num_iter = 200, report_step = 10)
  bsamp <- beta_auto[[1]]$sample_beta
  expect_equal(dim(bsamp), c(ncol(corr), 20))
  h2_est <- apply(bsamp, 2, function(x) crossprod(x, bigsparser::sp_prodVec(corr, x)))
  expect_equal(h2_est, beta_auto[[1]]$path_h2_est[seq(210, 400, by = 10)])

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
    snp_ldpred2_grid(corr, df_beta, params[c(1, 8), ], ncores = 2))
  all_beta <- replicate(
    n = 5, simplify = FALSE,
    snp_ldpred2_auto(corr, df_beta, h2_init = ldsc[["h2"]],
                     vec_p_init = seq_log(1e-4, 0.9, 6),
                     burn_in = 100, num_iter = 100,
                     sparse = TRUE, ncores = 2))
})

################################################################################
