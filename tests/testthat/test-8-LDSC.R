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
