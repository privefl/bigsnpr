################################################################################

context("SPLIT_LD")

################################################################################

test_that("get_L() and get_C() work", {

  corr <- as(outer(1:4 / 10, 1:4 / 10, '+'), "dgCMatrix")
  diag(corr) <- 1

  L <- matrix(0, 4, 4)
  L[1, 2] <- sum(corr[1, 2:4]^2)
  L[1, 3] <- sum(corr[1, 3:4]^2)
  L[1, 4] <- sum(corr[1, 4:4]^2)
  L[2, 3] <- sum(corr[2, 3:4]^2)
  L[2, 4] <- sum(corr[2, 4:4]^2)
  L[3, 4] <- sum(corr[3, 4:4]^2)

  # Rcpp::sourceCpp('src/split-LD.cpp')
  library(magrittr)
  L2 <- corr %>%
    Matrix::tril() %>%
    { get_L(.@p, .@i, .@x, thr_r2 = 0) } %>%
    # L now has an extra column with all 0s for convenience
    { Matrix::sparseMatrix(i = .[[1]], j = .[[2]], x = .[[3]], dims = c(4, 5),
                           triangular = FALSE, index1 = FALSE) }

  expect_equal(L2[, 1:4], as(L, "dgCMatrix"))

  # E <- matrix(0, 4, 4)
  # E[1, 1] <- sum(corr[  1, 2:4]^2)
  # E[1, 2] <- sum(corr[1:2, 3:4]^2)
  # E[1, 3] <- sum(corr[1:3,   4]^2)
  # E[2, 2] <- sum(corr[  2, 3:4]^2)
  # E[2, 3] <- sum(corr[2:3,   4]^2)
  # E[3, 3] <- sum(corr[  3,   4]^2)

  path_cost1 <- get_C(L2, min_size = 1, max_size = 4, K = 5)
  # one block remaining
  expect_identical(path_cost1$best_ind[, 1], rep(4L, 4))
  expect_identical(path_cost1$C[, 1], rep(0, 4))
  # two blocks remaining -> last is impossible
  expect_identical(path_cost1$best_ind[, 2], c(1L, 2L, 3L, NA))
  expect_equal(path_cost1$C[, 2], c(0.5, 0.61, 0.49, NA))
  # three blocks remaining
  expect_identical(path_cost1$best_ind[, 3], c(1L, 2L, NA, NA))
  expect_equal(path_cost1$C[, 3], c(1.11, 1.1, NA, NA), tolerance = 1e-6)
  # four blocks remaining
  expect_identical(path_cost1$best_ind[, 4], c(1L, NA, NA, NA))
  expect_equal(path_cost1$C[, 4], c(1.6, NA, NA, NA))
  # five blocks remaining -> all impossible
  expect_identical(path_cost1$best_ind[, 5], rep(NA_integer_, 4))
  expect_identical(path_cost1$C[, 5], rep(NA_real_, 4))

  # must use two blocks of two
  path_cost2 <- get_C(L2, min_size = 2, max_size = 2, K = 3)
  # one block remaining
  expect_identical(path_cost2$best_ind[, 1], c(NA, NA, 4L, NA))
  expect_identical(path_cost2$C[, 1], c(NA, NA, 0, NA))
  # two blocks remaining
  expect_identical(path_cost2$best_ind[, 2], c(2L, NA, NA, NA))
  expect_equal(path_cost2$C[, 2], c(1.02, NA, NA, NA), tolerance = 1e-6)
  # three blocks remaining -> always impossible
  expect_identical(path_cost2$best_ind[, 3], rep(NA_integer_, 4))
  expect_identical(path_cost2$C[, 3], rep(NA_real_, 4))

  path_cost3 <- get_C(L2, min_size = 1, max_size = 3, K = 3)
  # one block remaining
  expect_identical(path_cost3$best_ind[, 1], c(NA, 4L, 4L, 4L))
  expect_identical(path_cost3$C[, 1], c(NA, 0, 0, 0))
  # two blocks remaining
  expect_identical(path_cost3$best_ind[, 2], c(1L, 2L, 3L, NA))
  expect_equal(path_cost3$C[, 2], c(0.5, 0.61, 0.49, NA))
  # three blocks remaining
  expect_identical(path_cost3$best_ind[, 3], c(1L, 2L, NA, NA))
  expect_equal(path_cost3$C[, 3], c(1.11, 1.10, NA, NA), tolerance = 1e-6)
})

################################################################################

test_that("$perc_kept in snp_ldsplit() is exact", {

  corr <- as(outer(1:4 / 10, 1:4 / 10, '+'), "dgCMatrix")
  diag(corr) <- 1

  test <- snp_ldsplit(corr, 0, 1, 2, 4)
  expect_identical(test$n_block, 2:4)
  expect_identical(test$perc_kept, c(8, 6, 4) / 16)
})

################################################################################

test_that("snp_ldsplit() gives consistent results", {

  # skip_if(is_cran)

  corr <- readRDS(test_path("testdata/spMat.rds"))

  res1 <- snp_ldsplit(corr, thr_r2 = runif(1, 0, 0.2),
                      min_size = 10, max_size = 30, max_K = 50)
  # 30 * 13 = 390 < 401
  # 10 * 40 = 400 < 401
  expect_identical(res1$n_block, 14:40)

  res2 <- snp_ldsplit(corr, thr_r2 = runif(1, 0, 0.2),
                      min_size = 20, max_size = 40, max_K = 50)
  expect_identical(res2$n_block, 11:20)

  res3 <- snp_ldsplit(corr, thr_r2 = runif(1, 0, 0.2),
                      min_size = 20, max_size = 40, max_K = 15)
  expect_identical(res3$n_block, 11:15)
})

################################################################################
