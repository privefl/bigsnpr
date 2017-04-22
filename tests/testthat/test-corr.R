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

corr <- bigsnpr:::corMat(BM = attach.BM(G),
                         rowInd = ind.row,
                         colInd = ind.col,
                         size = size,
                         thr = rep(0.2, n))
corr2 <- corr^2

################################################################################

file.ld <- system.file("extdata", "example.ld", package = "bigsnpr")
true <- data.table::fread(file.ld, data.table = FALSE)
snps.ind <- sapply(true[, c("SNP_A", "SNP_B")], match,
                   table = paste0("SNP", ind.col - 1))

ind.size <- which(abs(snps.ind[, 1] - snps.ind[, 2]) <= size)

library(Matrix)
corr.true <- as(matrix(0, m, m), "dgCMatrix")
corr.true[snps.ind[ind.size, ]] <- true$R2[ind.size]

test_that("Same correlations as PLINK", {
  expect_equal(corr2@i,   corr.true@i)
  expect_equal(corr2@p,   corr.true@p)
  expect_equal(corr2@Dim, corr.true@Dim)
  expect_equal(corr2@x,   corr.true@x, tolerance = 1e-6)
})

################################################################################
