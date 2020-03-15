################################################################################

context("OPENMP")
# Basically, test if any crash..

skip_if(is_cran)

################################################################################

test_that("parallel snp_colstats() works", {

  G <- snp_attachExtdata()$genotypes
  rows <- sample(nrow(G), replace = TRUE)
  cols <- sample(ncol(G), replace = TRUE)

  test <- replicate(20, simplify = FALSE, {
    bigsnpr:::snp_colstats(G, rows, cols, ncores = 2)
  })
  true <- bigsnpr:::snp_colstats(G, rows, cols, ncores = 1)

  expect_true(all(sapply(test, all.equal, current = true)))
})

################################################################################

test_that("parallel bed_prodVec() works", {

  bedfile <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  obj.bed <- bed(bedfile)
  rows <- sample(nrow(obj.bed), replace = TRUE)
  cols <- sample(ncol(obj.bed), replace = TRUE)
  center <- rnorm(ncol(obj.bed))
  scale  <- runif(ncol(obj.bed))
  y.col <- rnorm(ncol(obj.bed))

  test <- replicate(20, simplify = FALSE, {
    bed_prodVec(obj.bed, y.col, rows, cols, center, scale, ncores = 2)
  })
  true <- bed_prodVec(obj.bed, y.col, rows, cols, center, scale, ncores = 1)

  expect_true(all(sapply(test, all.equal, current = true)))
})

################################################################################

test_that("parallel bed_cprodVec() works", {

  bedfile <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  obj.bed <- bed(bedfile)
  rows <- sample(nrow(obj.bed), replace = TRUE)
  cols <- sample(ncol(obj.bed), replace = TRUE)
  center <- rnorm(ncol(obj.bed))
  scale  <- runif(ncol(obj.bed))
  y.row <- rnorm(nrow(obj.bed))

  test <- replicate(20, simplify = FALSE, {
    bed_cprodVec(obj.bed, y.row, rows, cols, center, scale, ncores = 2)
  })
  true <- bed_cprodVec(obj.bed, y.row, rows, cols, center, scale, ncores = 1)

  expect_true(all(sapply(test, all.equal, current = true)))
})

################################################################################

test_that("parallel multLinReg()-bed works", {

  bedfile <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  obj.bed <- bed(bedfile)
  rows <- rows_along(obj.bed)
  cols <- sample(ncol(obj.bed), replace = TRUE)

  U <- bed_randomSVD(obj.bed, k = 3)$u

  test <- replicate(20, simplify = FALSE, {
    bigsnpr:::multLinReg(obj.bed, rows, cols, U, ncores = 2)
  })
  true <- bigsnpr:::multLinReg(obj.bed, rows, cols, U, ncores = 1)

  expect_true(all(sapply(test, all.equal, current = true)))
})

################################################################################

test_that("parallel multLinReg()-bigSNP works", {

  G <- snp_attachExtdata()$genotypes
  rows <- rows_along(G)
  cols <- sample(ncol(G), replace = TRUE)

  U <- big_SVD(G, fun.scaling = snp_scaleBinom(), k = 3)$u

  test <- replicate(20, simplify = FALSE, {
    bigsnpr:::multLinReg(G, rows, cols, U, ncores = 2)
  })
  true <- bigsnpr:::multLinReg(G, rows, cols, U, ncores = 1)

  expect_true(all(sapply(test, all.equal, current = true)))
})

################################################################################

test_that("parallel snp_clumping() works", {

  G <- snp_attachExtdata()$genotypes
  rows <- sample(nrow(G), replace = TRUE)

  test <- replicate(10, simplify = FALSE, {
    snp_clumping(G, infos.chr = rep(1, ncol(G)), ind.row = rows, ncores = 2)
  })
  true <- snp_clumping(G, infos.chr = rep(1, ncol(G)), ind.row = rows, ncores = 1)

  expect_true(all(sapply(test, identical, y = true)))
})

################################################################################

test_that("parallel snp_grid_clumping() works", {

  G <- snp_attachExtdata()$genotypes
  rows <- sample(nrow(G), replace = TRUE)
  lpS <- runif(ncol(G))

  test <- replicate(5, simplify = FALSE, {
    snp_grid_clumping(G, ind.row = rows, lpS = lpS, ncores = 2,
                      infos.chr = rep(1, ncol(G)), infos.pos = cols_along(G),
                      grid.base.size = 0.1, grid.thr.r2 = c(0.05, 0.5))
  })
  true <- snp_grid_clumping(G, ind.row = rows, lpS = lpS, ncores = 1,
                            infos.chr = rep(1, ncol(G)), infos.pos = cols_along(G),
                            grid.base.size = 0.1, grid.thr.r2 = c(0.05, 0.5))

  expect_true(all(sapply(test, identical, y = true)))
})

################################################################################

test_that("parallel snp_readBed2() works", {

  bedfile <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  obj.bed <- bed(bedfile)
  rows <- sample(nrow(obj.bed), replace = TRUE)
  cols <- sample(ncol(obj.bed), replace = TRUE)

  test <- replicate(20, simplify = FALSE, {
    snp_readBed2(bedfile, tempfile(), rows, cols, ncores = 2)
  })
  true <- snp_readBed2(bedfile, tempfile(), rows, cols, ncores = 2)

  expect_true(all(sapply(test, function(rds, current) {
    all.equal(snp_attach(rds)$genotypes[], current)
  }, current = snp_attach(true)$genotypes[])))
})

################################################################################

test_that("parallel snp_cor() works", {

  G <- snp_attachExtdata()$genotypes
  rows <- sample(nrow(G), replace = TRUE)
  cols <- sample(ncol(G), replace = TRUE)

  test <- replicate(5, simplify = FALSE, {
    snp_cor(G, rows, cols, ncores = 2)
  })
  true <- snp_cor(G, rows, cols, ncores = 1)

  expect_true(all(sapply(test, all.equal, current = true)))
})

################################################################################

test_that("parallel bed_counts() works", {

  bedfile <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  obj.bed <- bed(bedfile)
  rows <- sample(nrow(obj.bed), replace = TRUE)
  cols <- sample(ncol(obj.bed), replace = TRUE)

  test <- replicate(20, simplify = FALSE, {
    bed_counts(obj.bed, rows, cols, ncores = 2)
  })
  true <- bed_counts(obj.bed, rows, cols, ncores = 1)

  expect_true(all(sapply(test, identical, y = true)))
})

################################################################################

test_that("parallel snp_fastImputeSimple() works", {

  G <- snp_attachExtdata()$genotypes

  test <- replicate(20, simplify = FALSE, {
    GNA <- big_copy(G)
    ind <- sort(sample(length(GNA), length(GNA) / 100)); GNA[ind] <- 3
    method <- sample(c("mode", "mean0", "mean2", "random"), 1)
    G2 <- snp_fastImputeSimple(G, method = method, ncores = 2)
    sum(big_counts(G2)[4, ])
  })

  expect_true(all(test == 0))
})

################################################################################
