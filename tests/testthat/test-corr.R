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
t <- r * sqrt((length(ind.row)-2)/(1-r^2))

corr <- bigsnpr:::corMat(BM = G,
                         rowInd = ind.row,
                         colInd = ind.col,
                         size = size,
                         thr = rep(t, n))
corr2 <- corr^2

################################################################################

file.ld <- system.file("testdata", "example.ld", package = "bigsnpr")
true <- data.table::fread(file.ld, data.table = FALSE)
snps.ind <- sapply(true[, c("SNP_A", "SNP_B")], match,
                   table = paste0("SNP", ind.col - 1))
stopifnot(all(snps.ind[, 1] < snps.ind[, 2]))

ind.size <- which((snps.ind[, 2] - snps.ind[, 1]) <= size)

library(Matrix)
corr.true <- as(matrix(0, m, m), "dgCMatrix")
corr.true[snps.ind[ind.size, , drop = FALSE]] <- true$R2[ind.size]

test_that("Same correlations as PLINK", {
  expect_equal(corr2@i,   corr.true@i)
  expect_equal(corr2@p,   corr.true@p)
  expect_equal(corr2@Dim, corr.true@Dim)
  expect_equal(corr2@x,   corr.true@x, tolerance = 1e-6)
})

################################################################################

# Correlations with significance levels
suppressMessages(library(Hmisc))

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

true <- Hmisc::rcorr(G[ind.row, ind.col])
ind <- which(true$P < alpha)

corr <- snp_cor(G = G,
                ind.row = ind.row,
                ind.col = ind.col,
                alpha = alpha)

test_that("Same correlations as Hmisc", {
  expect_equal(dim(corr), c(m, m))
  expect_equal(corr[ind], true$r[ind])
  expect_equal(2*(length(corr@i) - nrow(corr)), length(ind))
})

# without diagonal
corr2 <- snp_cor(G = test$genotypes,
                 ind.row = ind.row,
                 ind.col = ind.col,
                 alpha = alpha,
                 fill.diag = FALSE)

test_that("Same correlations (no diagonal) as Hmisc", {
  expect_equal(dim(corr2), c(m, m))
  expect_equal(corr2[ind], true$r[ind])
  expect_equal(2*length(corr2@i), length(ind))
})

################################################################################
