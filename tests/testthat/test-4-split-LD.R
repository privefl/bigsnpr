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

  # can use one block
  path_cost1 <- get_C(L2, min_size = 1, max_size = 4, lambda = 0)
  expect_identical(path_cost1, list(C = rep(0, 5), best_ind = rep(4L, 4)))

  # must use two blocks of two
  path_cost2 <- get_C(L2, min_size = 2, max_size = 2, lambda = 0)
  expect_equal(path_cost2, list(C = c(1.02, NA, 0, NA, 0),
                                best_ind = c(2L, NA, 4L, NA)),
               tolerance = 1e-6)

  # have the choice between 1+3 (cost: 0.5) or 3+1 (cost: 1.1)
  path_cost3 <- get_C(L2, min_size = 1, max_size = 3, lambda = 0)
  expect_equal(path_cost3, list(C = c(0.5, 0, 0, 0, 0),
                                best_ind = c(1L, 4L, 4L, 4L)))
})

################################################################################

test_that("snp_ldsplit() gives consistent results", {

  skip_on_cran()
  skip_if_offline("dropbox.com")

  corr <- readRDS(url("https://www.dropbox.com/s/65u96jf7y32j2mj/spMat.rds?raw=1"))

  # Recall: results are now sorted by cost
  res1 <- snp_ldsplit(corr, 0, data.frame(min_size = c(5, 10, 20),
                                          max_size = 30,
                                          lambda = 0))
  expect_false(is.unsorted(res1$min_size))

  res2 <- snp_ldsplit(corr, 0, data.frame(min_size = 5,
                                          max_size = c(10, 20, 30),
                                          lambda = 0))
  expect_false(is.unsorted(res2$cost))
  expect_false(is.unsorted(rev(res2$max_size)))

  res3 <- snp_ldsplit(corr, 0, data.frame(min_size = 10,
                                          max_size = 20,
                                          lambda = c(0, 0.001, 0.01, 0.1)))
  expect_false(is.unsorted(res3$lambda))

  res4 <- snp_ldsplit(corr, 0.005, data.frame(min_size = 20,
                                              max_size = 30,
                                              lambda = 0.001))
  expect_length(res4$all_cost[[1]], ncol(corr) + 1L)
  path.end <- tail(res4$all_cost[[1]], 41)
  # starting from the end:
  #  - one 0 because C[m + 1]
  #  - first 19 are not possible because < min_size
  #  - from 20 to 30, it can be only one block (cannot be split)
  #  - from 31 to 39, they are not possible because > max_size
  #    and split would be < min_size
  expect_equal(rev(path.end[-1]), c(0, rep(c(NA, 0, NA), c(19, 11, 9))))
  expect_gt(path.end[1], 0)
})

################################################################################
