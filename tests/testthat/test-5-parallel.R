################################################################################

context("PARALLEL")

options(bigstatsr.check.parallel.blas = FALSE)

################################################################################

test_that("Sequential and Parallel", {

  test <- snp_attachExtdata()
  G <- test$genotypes
  CHR <- sort(sample(1:2, size = length(test$map$chromosome), replace = TRUE))
  POS <- test$map$physical.pos
  expect_equal(order(POS), seq_along(POS))

  expect_equal(snp_clumping(G, CHR),
               snp_clumping(G, CHR, ncores = 2))

  expect_equal(snp_autoSVD(G, CHR, POS, verbose = FALSE),
               snp_autoSVD(G, CHR, POS, ncores = NCORES, verbose = FALSE),
               tolerance = 1e-6)
})

################################################################################
