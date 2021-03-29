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
  seq_s <- sample(1:5 / 5, 4)
  beta_grid <- snp_lassosum2(corr, df_beta, s = seq_s,
                             nlambda = nlam, ncores = 2)
  expect_gt(max(cor(beta_grid, true_beta)), 0.4)
  params <- attr(beta_grid, "grid_param")
  expect_equal(nrow(params), 4 * nlam)
  expect_true(all(params$num_iter <= 501))
  expect_true(all(params$s %in% seq_s))

  # Errors
  expect_error(snp_lassosum2(corr, df_beta[-1]),
               "'df_beta' should have element 'beta'.")
})

################################################################################
