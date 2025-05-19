################################################################################

context("AUTO_SVD")

################################################################################

test <- snp_attachExtdata()
G   <- test$genotypes
CHR <- test$map$chromosome
POS <- test$map$physical.pos / 10
POS2 <- round(POS + 1)

obj.bed <- bed(snp_writeBed(test, bedfile = tempfile(fileext = ".bed")))

################################################################################

test_that("snp_autoSVD() works", {

  expect_error(snp_autoSVD(G, as.character(CHR), POS), "only integers")
  expect_error(snp_autoSVD(G, CHR[-1], POS2), bigstatsr:::GET_ERROR_DIM())
  expect_error(snp_autoSVD(G, CHR, POS2[-1]), bigstatsr:::GET_ERROR_DIM())
  expect_error(snp_autoSVD(G, CHR, min.mac = 0), "no variation; set min.mac > 0")
  expect_output(snp_autoSVD(G, CHR, POS2, thr.r2 = NA), "Skipping clumping.")

  expect_is(obj.svd <- snp_autoSVD(G, CHR, verbose = FALSE), "big_SVD")
  expect_identical(attr(obj.svd, "lrldr"), cbind(LD.wiki34[0, 1:3], Iter = integer()))

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

  obj.svd8 <- snp_autoSVD(G, CHR, POS, alpha.tukey = 0.9999, roll.size = 0,
                          int.min.size = 0, verbose = FALSE)
  lrldr8 <- attr(obj.svd8, "lrldr")
  expect_identical(sapply(lrldr8, typeof),
                   c(Chr = "integer", Start = "double", Stop = "double", Iter = "integer"))
  expect_gte(min(lrldr8$Iter), 1)

  obj.svd9 <- snp_autoSVD(G, CHR + 0, test$map$physical.pos, size = 5000, alpha.tukey = 0.9999,
                          roll.size = 0, int.min.size = 0, verbose = FALSE)
  lrldr9 <- attr(obj.svd9, "lrldr")
  expect_identical(sapply(lrldr9, typeof),
                   c(Chr = "double", Start = "integer", Stop = "integer", Iter = "integer"))
  expect_gte(min(lrldr9$Iter), 1)

  expect_output(
    snp_autoSVD(G, CHR, alpha.tukey = 0.999999999, roll.size = 0, verbose = TRUE),
    "Maximum number of iterations reached.")
})

################################################################################

test_that("bed_autoSVD() works", {

  expect_error(bed_autoSVD(obj.bed, min.mac = 0), "no variation; set min.mac > 0")
  expect_output(bed_autoSVD(obj.bed, thr.r2 = NA), "Skipping clumping.")

  expect_is(obj.svd <- bed_autoSVD(obj.bed, verbose = FALSE), "big_SVD")
  expect_identical(attr(obj.svd, "lrldr"), cbind(LD.wiki34[0, 1:3], Iter = integer()))

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

test_that("MAC/MAF thresholds work", {

  info <- bed_MAF(obj.bed)

  replicate(10, {

    min.mac <- sample(1:40, size = 1)
    min.maf <- runif(n = 1, min = 0.01, max = 0.1)

    obj.svd1 <- snp_autoSVD(G, CHR, size = 5, min.mac = min.mac, min.maf = min.maf,
                            thr.r2 = NA, max.iter = 0, verbose = FALSE)
    ind1 <- attr(obj.svd1, "subset")
    expect_true(all(info$maf[ind1] >= min.maf))
    expect_true(all(info$mac[ind1] >= min.mac))

    obj.svd2 <- bed_autoSVD(obj.bed, min.mac = min.mac, min.maf = min.maf,
                            thr.r2 = NA, max.iter = 0, verbose = FALSE)
    ind2 <- attr(obj.svd2, "subset")
    expect_true(all(info$maf[ind2] >= min.maf))
    expect_true(all(info$mac[ind2] >= min.mac))
  })
})

################################################################################
