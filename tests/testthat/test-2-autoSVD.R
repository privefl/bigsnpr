################################################################################

context("AUTO_SVD")

################################################################################

test <- snp_attachExtdata()
G   <- test$genotypes
CHR <- test$map$chromosome
POS <- test$map$physical.pos / 10
POS2 <- round(POS + 1)

################################################################################

test_that("snp_autoSVD() works", {

  expect_error(snp_autoSVD(G, CHR, POS), "only integers")
  expect_error(snp_autoSVD(G, CHR[-1], POS2), bigstatsr:::GET_ERROR_DIM())
  expect_error(snp_autoSVD(G, CHR, POS2[-1]), bigstatsr:::GET_ERROR_DIM())
  expect_output(snp_autoSVD(G, CHR, POS2, thr.r2 = NA), "Skipping clumping.")

  expect_is(obj.svd <- snp_autoSVD(G, CHR, verbose = FALSE), "big_SVD")
  expect_identical(attr(obj.svd, "lrldr"), LD.wiki34[0, 1:3])

  skip_if(is_cran)

  obj.svd2 <- snp_autoSVD(G, CHR, size = 5, verbose = FALSE)
  expect_gt(length(attr(obj.svd2, "subset")), length(attr(obj.svd, "subset")))
  obj.svd3 <- snp_autoSVD(G, CHR, POS2, size = 5, verbose = FALSE)
  expect_lt(length(attr(obj.svd3, "subset")), length(attr(obj.svd2, "subset")))

  obj.svd4 <- snp_autoSVD(G, CHR, roll.size = 0, verbose = FALSE)
  expect_lt(length(attr(obj.svd4, "subset")), length(attr(obj.svd, "subset")))

  obj.svd5 <- snp_autoSVD(G, CHR, thr.r2 = 1, roll.size = 0, verbose = FALSE)
  expect_gt(min(diag(cor(obj.svd5$u, obj.svd$u))[1:3]^2), 0.98)

  obj.svd6 <- snp_autoSVD(G, CHR, thr.r2 = NA, roll.size = 0, verbose = FALSE)
  expect_equal(obj.svd6, obj.svd5)

  obj.svd7 <- snp_autoSVD(G, CHR, alpha.tukey = 0.999, roll.size = 0, verbose = FALSE)
  expect_lt(length(attr(obj.svd7, "subset")), length(attr(obj.svd6, "subset")))

  expect_output(
    snp_autoSVD(G, CHR, alpha.tukey = 0.999999999, roll.size = 0, verbose = TRUE),
    "Maximum number of iterations reached.")
})

################################################################################

test_that("bed_autoSVD() works", {

  obj.bed <- bed(snp_writeBed(test, bedfile = tempfile(fileext = ".bed")))

  expect_output(bed_autoSVD(obj.bed, thr.r2 = NA), "Skipping clumping.")

  expect_is(obj.svd <- bed_autoSVD(obj.bed, verbose = FALSE), "big_SVD")
  expect_identical(attr(obj.svd, "lrldr"), LD.wiki34[0, 1:3])

  skip_if(is_cran)

  obj.svd2 <- bed_autoSVD(obj.bed, size = 5, verbose = FALSE)
  expect_gt(length(attr(obj.svd2, "subset")), length(attr(obj.svd, "subset")))

  obj.svd4 <- bed_autoSVD(obj.bed, roll.size = 0, verbose = FALSE)
  expect_lt(length(attr(obj.svd4, "subset")), length(attr(obj.svd, "subset")))

  obj.svd5 <- bed_autoSVD(obj.bed, thr.r2 = 1, roll.size = 0, verbose = FALSE)
  expect_gt(min(diag(cor(obj.svd5$u, obj.svd$u))[1:3]^2), 0.98)

  obj.svd6 <- bed_autoSVD(obj.bed, thr.r2 = NA, roll.size = 0, verbose = FALSE)
  expect_equal(obj.svd6, obj.svd5)

  obj.svd7 <- bed_autoSVD(obj.bed, alpha.tukey = 0.999, roll.size = 0, verbose = FALSE)
  expect_lt(length(attr(obj.svd7, "subset")), length(attr(obj.svd6, "subset")))

  expect_output(
    bed_autoSVD(obj.bed, alpha.tukey = 0.999999999, roll.size = 0, verbose = TRUE),
    "Maximum number of iterations reached.")
})

################################################################################
