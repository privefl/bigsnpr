################################################################################

context("CORR")

################################################################################

test <- snp_attachExtdata()
G <- test$genotypes
n <- nrow(G)
m <- ncol(G)
ind.row <- seq_len(n)
ind.col <- seq_len(m)

size <- round(runif(1, 1, 200))

r <- sqrt(0.2)
t <- r * sqrt((length(ind.row) - 2) / (1 - r^2))
r_again <- t / sqrt(length(ind.row) - 2 + t^2)
expect_equal(r_again, r)

ind_val <- bigsnpr:::corMat(G,
                            rowInd = ind.row,
                            colInd = ind.col,
                            size = size,
                            thr = rep(r, n),
                            pos = seq_along(ind.col),
                            ncores = 2)

list_i <- lapply(ind_val, function(.) .$i)
m <- length(list_i)
seq_ind <- seq_len(m)
i <- unlist(list_i)
j <- rep(seq_ind, lengths(list_i))
x <- unlist(lapply(ind_val, function(.) .$x))

corr2 <- Matrix::sparseMatrix(i = i, j = j, x = x ** 2, dims = c(m, m))

################################################################################

file.ld <- system.file("testdata", "example.ld", package = "bigsnpr")
true <- bigreadr::fread2(file.ld)
snps.ind <- sapply(true[, c("SNP_A", "SNP_B")], match,
                   table = paste0("SNP", ind.col - 1))
stopifnot(all(snps.ind[, 1] < snps.ind[, 2]))

ind.size <- which((snps.ind[, 2] - snps.ind[, 1]) <= size)

library(Matrix)
corr.true <- Matrix(0, m, m, sparse = TRUE, doDiag = FALSE)
corr.true[snps.ind[ind.size, , drop = FALSE]] <- true$R2[ind.size]

test_that("Same correlations as PLINK", {
  expect_equal(corr2@i,   corr.true@i)
  expect_equal(corr2@p,   corr.true@p)
  expect_equal(corr2@Dim, corr.true@Dim)
  expect_equal(corr2@x,   corr.true@x, tolerance = 1e-6)
})

################################################################################

N <- 500
M <- 100
test <- snp_fake(N, M)
G <- test$genotypes
G[] <- sample(as.raw(0:3), size = length(G), replace = TRUE)
corr_diag <- snp_cor(G, thr_r2 = 0.99)
expect_identical(corr_diag@uplo, "U")  # case of diagonal matrix
expect_s4_class(as_SFBM(corr_diag), "SFBM")

# random parameters
alpha <- runif(1, 0.01, 0.2)
ind.row <- sample(N, N / 2)
ind.col <- sample(M, M / 2)
m <- length(ind.col)

test_that("Same correlations as Hmisc (with significance levels)", {

  skip_if_not_installed("Hmisc")

  true <- Hmisc::rcorr(G[ind.row, ind.col])
  ind <- which(true$P < alpha)

  corr <- snp_cor(Gna = G,
                  ind.row = ind.row,
                  ind.col = ind.col,
                  alpha = alpha)

  expect_equal(dim(corr), c(m, m))
  expect_equal(corr[ind], true$r[ind])
  expect_equal(2*(length(corr@i) - nrow(corr)), length(ind))

  # without diagonal
  corr2 <- snp_cor(G = G,
                   ind.row = ind.row,
                   ind.col = ind.col,
                   alpha = alpha,
                   fill.diag = FALSE)

  expect_equal(dim(corr2), c(m, m))
  expect_equal(corr2[ind], true$r[ind])
  expect_equal(2*length(corr2@i), length(ind))

  # with additional threshold
  corr3 <- snp_cor(G = G,
                   ind.row = ind.row,
                   ind.col = ind.col,
                   alpha = alpha,
                   thr_r2 = 0.02,
                   fill.diag = TRUE)

  expect_equal(dim(corr3), c(m, m))
  expect_gte(min(corr3@x^2), 0.02)
})

################################################################################

test_that("Information on position is used in snp_cor()", {

  corr3 <- snp_cor(G = G, ind.row = ind.row, ind.col = ind.col,
                   alpha = alpha, fill.diag = FALSE, size = 5)

  corr3.2 <- snp_cor(G = G, ind.row = ind.row, ind.col = ind.col,
                    alpha = alpha, fill.diag = FALSE, size = 5,
                    ncores = 2)
  expect_equal(corr3.2, corr3)

  corr4 <- snp_cor(G = G, ind.row = ind.row, ind.col = ind.col,
                   alpha = alpha, fill.diag = FALSE, size = 5e3,
                   infos.pos = 1e6 * seq_along(ind.col))
  expect_equal(corr4, corr3)

  corr5 <- snp_cor(G = G, ind.row = ind.row, ind.col = ind.col,
                   alpha = alpha, fill.diag = FALSE, size = 5e-3,
                   infos.pos = seq_along(ind.col))
  expect_equal(corr5, corr3)

  corr6 <- snp_cor(G = G, ind.row = ind.row, ind.col = ind.col,
                   alpha = alpha, fill.diag = FALSE, size = 5e-3,
                   infos.pos = 1000 * seq_along(ind.col))
  expect_length(corr6@x, 0)
})

################################################################################

test_that("bed_cor() works like snp_cor()", {

  bedfile <- snp_writeBed(test, tempfile(fileext = ".bed"))
  obj_bed <- bed(bedfile)
  alpha <- runif(1)
  expect_error(bed_cor(obj.bed = G, alpha = alpha),
               "'obj.bed' is not of class 'bed'.", fixed = TRUE)
  expect_equal(
    snp_cor(Gna = G, alpha = alpha),
    bed_cor(obj.bed = obj_bed, alpha = alpha)
  )
})

################################################################################

test_that("snp_cor() returns NaNs when SD is 0 (and warns about it)", {

  G <- snp_fake(10, 10)$genotypes
  G[] <- sample(0:2, length(G), replace = TRUE)
  G[, 1] <- 0
  corr <- expect_warning(snp_cor(G), "NA or NaN values")
  corr2 <- expect_warning(cor(G[]), "the standard deviation is zero")
  expect_equal(as.matrix(corr), corr2, check.attributes = FALSE)
})

################################################################################
