################################################################################

context("LASSOSUM2")

################################################################################

test_that("lassosum2 works", {

  skip_if(is_cran)
  skip_if_offline("dropbox.com")

  load(url("https://www.dropbox.com/s/c13uygnjh6yh7vf/to-test-ldpred2.RData?raw=1"))
  corr <- as_SFBM(corr)

  # lassosum2
  nlam <- sample(8:15, 1)
  beta_grid <- snp_lassosum2(corr, df_beta, nlambda = nlam, ncores = 2)
  expect_gt(max(cor(beta_grid, true_beta), na.rm = TRUE), 0.4)
  params <- attr(beta_grid, "grid_param")
  expect_equal(nrow(params), 6 * nlam)
  expect_true(all(params$num_iter <= 501))

  # Errors
  expect_error(snp_lassosum2(corr, df_beta[-1]),
               "'df_beta' should have element 'beta'.")
  expect_error(snp_lassosum2(corr, df_beta[-2]),
               "'df_beta' should have element 'beta_se'.")
})

################################################################################
