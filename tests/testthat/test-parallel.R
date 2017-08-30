################################################################################

context("PARALLEL")

# https://github.com/hadley/testthat/issues/567
Sys.unsetenv("R_TESTS")

################################################################################

test <- snp_attachExtdata()

G <- test$genotypes
CHR <- sort(sample(1:2, size = length(test$map$chromosome), replace = TRUE))
POS <- test$map$physical.pos
expect_equal(order(POS), seq_along(POS))

################################################################################

test_that("Sequential and Parallel", {
  expect_equal(snp_pruning(G, CHR),
               snp_pruning(G, CHR, ncores = 2))

  expect_equal(snp_clumping(G, CHR),
               snp_clumping(G, CHR, ncores = 2))

  skip_on_cran()
  expect_equal(snp_autoSVD(G, CHR, POS),
               snp_autoSVD(G, CHR, POS, ncores = 2), tolerance = 1e-6)
})

################################################################################

test_that("Sequential and Parallel + seeds", {
  elemNA <- sample(length(G), size = 20)

  G2 <- big_copy(G); G2[elemNA] <- as.raw(3)  # NA
  G3 <- big_copy(G); G3[elemNA] <- as.raw(3)  # NA

  expect_equal(snp_fastImpute(G2, CHR, seed = 2),
               snp_fastImpute(G3, CHR, seed = 2, ncores = 2))

  expect_equal(G2[], G3[])
})

################################################################################
