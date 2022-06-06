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
    { bigsnpr:::get_L(.@p, .@i, .@x, thr_r2 = 0, max_r2 = 1) } %>%
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

  path_cost1 <- bigsnpr:::get_C(L2, min_size = 1, max_size = 4,
                                max_K = 5, max_cost = Inf)
  # one block remaining
  expect_identical(path_cost1$best_ind[, 1], rep(4L, 4))
  expect_identical(path_cost1$C[, 1], rep(0, 4))
  # two blocks remaining -> last is impossible
  expect_identical(path_cost1$best_ind[, 2], c(1L, 2L, 3L, NA))
  expect_equal(path_cost1$C[, 2], c(0.5, 0.61, 0.49, Inf))
  # three blocks remaining
  expect_identical(path_cost1$best_ind[, 3], c(1L, 2L, NA, NA))
  expect_equal(path_cost1$C[, 3], c(1.11, 1.1, Inf, Inf), tolerance = 1e-6)
  # four blocks remaining
  expect_identical(path_cost1$best_ind[, 4], c(1L, NA, NA, NA))
  expect_equal(path_cost1$C[, 4], c(1.6, Inf, Inf, Inf))
  # five blocks remaining -> all impossible
  expect_identical(path_cost1$best_ind[, 5], rep(NA_integer_, 4))
  expect_identical(path_cost1$C[, 5], rep(Inf, 4))

  # must use two blocks of two
  path_cost2 <- bigsnpr:::get_C(L2, min_size = 2, max_size = 2,
                                max_K = 3, max_cost = Inf)
  # one block remaining
  expect_identical(path_cost2$best_ind[, 1], c(NA, NA, 4L, NA))
  expect_identical(path_cost2$C[, 1], c(Inf, Inf, 0, Inf))
  # two blocks remaining
  expect_identical(path_cost2$best_ind[, 2], c(2L, NA, NA, NA))
  expect_equal(path_cost2$C[, 2], c(1.02, Inf, Inf, Inf), tolerance = 1e-6)
  # three blocks remaining -> always impossible
  expect_identical(path_cost2$best_ind[, 3], rep(NA_integer_, 4))
  expect_identical(path_cost2$C[, 3], rep(Inf, 4))

  path_cost3 <- bigsnpr:::get_C(L2, min_size = 1, max_size = 3,
                                max_K = 3, max_cost = Inf)
  # one block remaining
  expect_identical(path_cost3$best_ind[, 1], c(NA, 4L, 4L, 4L))
  expect_identical(path_cost3$C[, 1], c(Inf, 0, 0, 0))
  # two blocks remaining
  expect_identical(path_cost3$best_ind[, 2], c(1L, 2L, 3L, NA))
  expect_equal(path_cost3$C[, 2], c(0.5, 0.61, 0.49, Inf))
  # three blocks remaining
  expect_identical(path_cost3$best_ind[, 3], c(1L, 2L, NA, NA))
  expect_equal(path_cost3$C[, 3], c(1.11, 1.10, Inf, Inf), tolerance = 1e-6)
})

################################################################################

test_that("$perc_kept in snp_ldsplit() is exact", {

  corr <- as(outer(1:4 / 10, 1:4 / 10, '+'), "dgCMatrix")
  diag(corr) <- 1

  test <- snp_ldsplit(corr, 0, 1, 2, 4, max_r2 = 1, max_cost = Inf)
  expect_identical(test$n_block, 2:4)
  expect_identical(test$perc_kept, c(8, 6, 4) / 16)
})

################################################################################

test_that("snp_ldsplit() gives consistent results", {

  compute_cost <- function(block_num, corr, thr_r2) {
    corr %>%
      Matrix::tril() %>%
      as("dgTMatrix") %>%
      { .@x[block_num[.@i + 1L] != block_num[.@j + 1L]]^2 } %>%
      { sum(.[. >= thr_r2]) }
  }

  compute_max_out <- function(block_num, corr) {
    corr %>%
      Matrix::tril() %>%
      as("dgTMatrix") %>%
      { .@x[block_num[.@i + 1L] != block_num[.@j + 1L]]^2 } %>%
      max()
  }

  corr <- readRDS(test_path("testdata/spMat.rds"))

  res1 <- snp_ldsplit(corr, thr_r2 = 0.02, min_size = 10, max_size = 30,
                      max_K = 50, max_r2 = 1, max_cost = Inf)
  # 30 * 13 = 390 < 401
  # 10 * 40 = 400 < 401
  expect_identical(res1$n_block, 14:40)

  block_num <- lapply(res1$all_size, function(sizes) rep(seq_along(sizes), sizes))
  costs <- sapply(block_num, function(nums) compute_cost(nums, corr, thr_r2 = 0.02))
  expect_equal(res1$cost, costs)

  prev_res1 <- readRDS(test_path("testdata/split_before.rds"))  # < v1.10.1
  expect_equal(res1$cost, prev_res1$cost)  # only costs are the same

  res2 <- snp_ldsplit(corr, thr_r2 = runif(1, 0, 0.2), min_size = 20, max_size = 40,
                      max_K = 50, max_r2 = 1, max_cost = Inf)
  expect_identical(res2$n_block, 11:20)

  res3 <- snp_ldsplit(corr, thr_r2 = runif(1, 0, 0.2), min_size = 20, max_size = 40,
                      max_K = 15, max_r2 = 1, max_cost = Inf)
  expect_identical(res3$n_block, 11:15)

  # parameter 'max_cost' works
  max_cost <- runif(1, min(res1$cost), max(res1$cost))
  res4 <- snp_ldsplit(corr, thr_r2 = 0.02, min_size = 10, max_size = 30,
                      max_K = 50, max_r2 = 1, max_cost = max_cost)
  expect_true(all(res4$cost <= max_cost))
  res1.bad <- subset(res1, !n_block %in% res4$n_block)
  expect_true(all(res1.bad$cost > max_cost))

  # parameter 'max_r2' works
  max_r2 <- runif(1, 0.1, 0.4)
  res5 <- snp_ldsplit(corr, thr_r2 = 0.02, min_size = 10, max_size = 50,
                      max_K = 100, max_r2 = max_r2, max_cost = Inf)
  sapply(res5$all_size, function(sizes) {
     block_num <- rep(seq_along(sizes), sizes)
     max_r2_out <- compute_max_out(block_num, corr)
     expect_lte(max_r2_out, max_r2)
     max_r2_out / max_r2
  })

  # can provide multiple 'max_size'
  max_r2 <- runif(1, 0.3, 1)
  expect_equal(
    rbind(snp_ldsplit(corr, thr_r2 = 0.02, min_size = 10, max_size = 30,
                      max_K = 50, max_r2 = max_r2, max_cost = Inf),
          snp_ldsplit(corr, thr_r2 = 0.02, min_size = 10, max_size = 40,
                      max_K = 50, max_r2 = max_r2, max_cost = Inf)),
    snp_ldsplit(corr, thr_r2 = 0.02, min_size = 10, max_size = c(30, 40),
                max_K = 50, max_r2 = max_r2, max_cost = Inf)
  )
})

################################################################################
