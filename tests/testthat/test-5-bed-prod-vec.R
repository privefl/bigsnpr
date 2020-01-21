################################################################################

context("BED_PROD_VEC")

read_bed_scaled <- bigsnpr:::read_bed_scaled
ERROR_DIM <- bigstatsr:::GET_ERROR_DIM()

################################################################################

bedfile <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
obj.bed <- bed(bedfile)

N <- nrow(obj.bed)
M <- ncol(obj.bed)

################################################################################

test_that("equality with %*%", {
  replicate(20, {
    n <- sample(N, size = 1)
    m <- sample(M, size = 1)
    ind.row <- sample(N, size = n)
    ind.col <- sample(M, size = m)
    center <- rep(0, m)
    scale <- rep(1, m)
    y.col <- rnorm(m)
    X <- read_bed_scaled(obj.bed, ind.row, ind.col, center, scale)
    expect_equal(bed_prodVec(obj.bed, y.col, ind.row, ind.col),
                 drop(X %*% y.col))
    y.row <- rnorm(n)
    expect_equal(bed_cprodVec(obj.bed, y.row, ind.row, ind.col),
                 drop(crossprod(X, y.row)))

    center <- rnorm(m); scale <- runif(m)
    X <- read_bed_scaled(obj.bed, ind.row, ind.col, center, scale)
    expect_equal(bed_prodVec(obj.bed, y.col, ind.row, ind.col, center, scale),
                 drop(X %*% y.col))
    expect_equal(bed_cprodVec(obj.bed, y.row, ind.row, ind.col, center, scale),
                 drop(crossprod(X, y.row)))
  })
})

test_that("Incompatiblity between dimensions", {
  ind.row <- sample(N, size = 21)
  ind.col <- sample(M, size = 11)
  y.col <- rnorm(21)
  expect_error(bed_prodVec(obj.svd, y.col, ind.row, ind.col), ERROR_DIM)
  y.row <- rnorm(11)
  expect_error(bed_cprodVec(obj.bed, y.row, ind.row, ind.col), ERROR_DIM)
})

################################################################################
