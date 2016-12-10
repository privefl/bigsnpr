################################################################################

context("PRUNE")

test <- snp_readExample("test-pruning")

################################################################################

test$fam$pheno <- sample(1:2, size = nrow(test$fam),
                         replace = TRUE)
ind <- snp_pruning(test)
ind3 <- snp_pruning(test, thr.corr = 0.2)

################################################################################

test_that("Lengths differ with thresholds", {
  expect_gte(length(ind), length(ind3))
})

################################################################################
