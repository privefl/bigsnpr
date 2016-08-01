################################################################################

context("AUC")

x <- rnorm(134)
y <- sample(c(-1, 1), size = length(x), replace = TRUE)

auc <- AucSample(x, y, nsim = 1e6, seed = 1)
auc2 <- AucSample(x, y, nsim = 1e6, seed = 1)
auc.conf <- AucSampleConf(x, y, nboot = 1e4, nsim = 1e3, seed = 1)
auc.conf2 <- AucSampleConf(x, y, nboot = 1e4, nsim = 1e3, seed = 1)

################################################################################

test_that("Same results of AUC with seed", {
  expect_equal(auc, auc2)
  expect_equal(auc.conf, auc.conf2)
})

################################################################################
