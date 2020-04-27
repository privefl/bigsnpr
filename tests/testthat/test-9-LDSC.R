################################################################################

context("LDSC")

################################################################################

test_that("wlm() works", {

  x <- rnorm(100)
  y <- rnorm(100)
  w <- runif(100)

  my_wlm <- stats::lm.wfit(cbind(1, x), y, w)
  my_wlm2 <- wlm(x, y, w)
  expect_equal(my_wlm2$intercept, my_wlm$coefficients[1], check.attributes = FALSE)
  expect_equal(my_wlm2$slope,     my_wlm$coefficients[2], check.attributes = FALSE)
  expect_equal(my_wlm2$pred,      my_wlm$fitted.values)

  my_wlm <- stats::lm.wfit(as.matrix(x), y, w)
  my_wlm2 <- wlm_no_int(x, y, w)
  expect_equal(my_wlm2$slope, my_wlm$coefficients, check.attributes = FALSE)
  expect_equal(my_wlm2$pred,  my_wlm$fitted.values)
})

################################################################################

test_that("snp_ldsc() gives the same result as bulik/ldsc", {

  skip_if(is_cran)
  skip_if_offline("dropbox.com")

  input <- readRDS(url(
    "https://www.dropbox.com/s/fhebnl9dn9poevp/sumstats-ldsc.rds?raw=1"))

  input$ncores <- 2
  ldsc <- do.call(bigsnpr::snp_ldsc, args = input)
  # 0.772006663 0.008230629 0.140788419 0.006605051
  expect_equal(ldsc, c(0.772, 0.0082, 0.1408, 0.0066),
               check.attributes = FALSE, tolerance = 4e-5)

  input$intercept <- 1
  ldsc_no_int <- do.call(bigsnpr::snp_ldsc, args = input)
  # 1.000000000 0.000000000 0.039032515 0.004779765
  expect_equal(ldsc_no_int[1:3], c(1, 0, 0.039),
               check.attributes = FALSE, tolerance = 4e-5)
  expect_equal(ldsc_no_int[[4]], 0.0057, tolerance = 1e-3)
  # not really the same -> due to the fast block jackknives they use?

  input$intercept <- NULL
  input$chi2_thr1 <- 10
  ldsc_max_chi2 <- do.call(bigsnpr::snp_ldsc, args = input)
  # 0.777193737 0.006735874 0.138051207 0.006193540
  expect_equal(ldsc_max_chi2, c(0.7772, 0.0069, 0.1381, 0.0061),
               check.attributes = FALSE, tolerance = 2e-4)
})

################################################################################
