################################################################################

context("SPLIT_LD")

################################################################################

test_that("get_L() and get_E() work", {

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
    { get_L(.@p, .@i, .@x) } %>%
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

  E2 <- L2 %>%
    get_E(min_row = min_row(corr@p, corr@i)) %>%
    { Matrix::sparseMatrix(i = .[[1]], j = .[[2]], x = .[[3]], dims = c(4, 4),
                           triangular = FALSE, index1 = TRUE) }

  expect_equal(E2, as(E, "dgCMatrix"))
})

################################################################################
