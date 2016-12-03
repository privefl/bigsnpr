################################################################################

context("AUC")

x <- rnorm(134)
y <- sample(c(-1, 1), size = length(x), replace = TRUE)

auc <- snp_aucSample(x, y, nsim = 1e6, seed = 1)
auc2 <- snp_aucSample(x, y, nsim = 1e6, seed = 1)
auc.conf <- snp_aucSampleConf(x, y, nboot = 1e4, nsim = 1e3, seed = 1)
auc.conf2 <- snp_aucSampleConf(x, y, nboot = 1e4, nsim = 1e3, seed = 1)

################################################################################

test_that("Same results of AUC with seed", {
  expect_equal(auc, auc2)
  expect_equal(auc.conf, auc.conf2)
})

################################################################################

test_that("Same results of AUC in particular cases", {
  expect_equal(snp_aucSample(c(0, 0), 0:1), 0.5) # Equality of scores
  expect_equal(snp_aucSample(c(0.2, 0.1, 1), c(-1, -1, 1)), 1) # Perfect AUC
  expect_equivalent(snp_aucSampleConf(c(0, 0), 0:1, nboot = 1e3),
                    c(rep(0.5, 3), 0))
  expect_equivalent(
    snp_aucSampleConf(c(0.2, 0.1, 1), c(-1, -1, 1), nboot = 1e3),
                    c(rep(1, 3), 0))
})

################################################################################
