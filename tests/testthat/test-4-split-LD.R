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
    { Matrix::sparseMatrix(i = .[[1]], j = .[[2]], x = .[[3]], dims = c(4, 4),
                           triangular = FALSE, index1 = FALSE) }

  expect_equal(L2, as(L, "dgCMatrix"))

  E <- matrix(0, 4, 4)
  E[1, 1] <- sum(corr[  1, 2:4]^2)
  E[1, 2] <- sum(corr[1:2, 3:4]^2)
  E[1, 3] <- sum(corr[1:3,   4]^2)
  E[2, 2] <- sum(corr[  2, 3:4]^2)
  E[2, 3] <- sum(corr[2:3,   4]^2)
  E[3, 3] <- sum(corr[  3,   4]^2)

  # can use one block
  path_cost1 <- get_C(L2, min_size = 1, max_size = 4, lambda = 0)
  expect_identical(path_cost1, list(C = rep(0, 4), best_ind = rep(NA_integer_, 4)))

  # must use two blocks of two
  path_cost2 <- get_C(L2, min_size = 2, max_size = 2, lambda = 0)
  expect_equal(path_cost2, list(C = c(1.02, NA, 0, NA),
                                best_ind = c(2L, NA, NA, NA)),
               tolerance = 1e-6)

  # have the choice between 1+3 (cost: 0.5) or 3+1 (cost: 1.1)
  path_cost3 <- get_C(L2, min_size = 1, max_size = 3, lambda = 0)
  expect_equal(path_cost3, list(C = c(0.5, 0, 0, 0),
                                best_ind = c(1L, NA, NA, NA)))
})

################################################################################

test_that("snp_ldsplit() gives consistent results", {

  skip_on_cran()
  skip_if_offline("dropbox.com")

  corr <- readRDS(url("https://www.dropbox.com/s/65u96jf7y32j2mj/spMat.rds?raw=1"))

  res1 <- snp_ldsplit(corr, 0, data.frame(min_size = c(5, 10, 20),
                                          max_size = 30,
                                          lambda = 0))
  expect_false(is.unsorted(res1$cost))


  res2 <- snp_ldsplit(corr, 0, data.frame(min_size = 5,
                                          max_size = c(10, 20, 30),
                                          lambda = 0))
  expect_false(is.unsorted(rev(res2$cost)))

  res3 <- snp_ldsplit(corr, 0, data.frame(min_size = 10,
                                          max_size = 20,
                                          lambda = c(0, 0.001, 0.1)))
  expect_false(is.unsorted(res3$cost))
})

################################################################################
