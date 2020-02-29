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

corr <- bigsnpr:::corMat(BM = G,
                         rowInd = ind.row,
                         colInd = ind.col,
                         blockInd = seq_along(ind.col),
                         size = size,
                         thr = rep(r, n),
                         pos = seq_along(ind.col))
corr2 <- corr^2

################################################################################

file.ld <- system.file("testdata", "example.ld", package = "bigsnpr")
true <- bigreadr::fread2(file.ld)
snps.ind <- sapply(true[, c("SNP_A", "SNP_B")], match,
                   table = paste0("SNP", ind.col - 1))
stopifnot(all(snps.ind[, 1] < snps.ind[, 2]))

ind.size <- which((snps.ind[, 2] - snps.ind[, 1]) <= size)

library(Matrix)
corr.true <- Matrix(0, m, m, sparse = TRUE)
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
  corr2 <- snp_cor(G = test$genotypes,
                   ind.row = ind.row,
                   ind.col = ind.col,
                   alpha = alpha,
                   fill.diag = FALSE)

  expect_equal(dim(corr2), c(m, m))
  expect_equal(corr2[ind], true$r[ind])
  expect_equal(2*length(corr2@i), length(ind))
})

################################################################################

test_that("Information on position is used in snp_cor()", {

  corr3 <- snp_cor(G = test$genotypes, ind.row = ind.row, ind.col = ind.col,
                   alpha = alpha, fill.diag = FALSE, size = 5)

  corr3.2 <- snp_cor(G = test$genotypes, ind.row = ind.row, ind.col = ind.col,
                    alpha = alpha, fill.diag = FALSE, size = 5,
                    ncores = 2)
  expect_equal(corr3.2, corr3)

  corr4 <- snp_cor(G = test$genotypes, ind.row = ind.row, ind.col = ind.col,
                   alpha = alpha, fill.diag = FALSE, size = 5e3,
                   infos.pos = 1e6 * seq_along(ind.col))
  expect_equal(corr4, corr3)

  corr5 <- snp_cor(G = test$genotypes, ind.row = ind.row, ind.col = ind.col,
                   alpha = alpha, fill.diag = FALSE, size = 5e-3,
                   infos.pos = seq_along(ind.col))
  expect_equal(corr5, corr3)

  corr6 <- snp_cor(G = test$genotypes, ind.row = ind.row, ind.col = ind.col,
                   alpha = alpha, fill.diag = FALSE, size = 5e-3,
                   infos.pos = 1000 * seq_along(ind.col))
  expect_length(corr6@x, 0)
})

################################################################################
