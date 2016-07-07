################################################################################

context("REGLIN")

opt.save <- options(bigmemory.typecast.warning = FALSE)

# Simulating some data
N1 <- 3000
N2 <- 1000
x1 <- rnorm(N1, 0.5)
x2 <- rnorm(N2, -0.5)
x <- c(x1, x2)

################################################################################

get_res_weights <- function(x.big, y) {
  res <- matrix(0, 3, ncol(x.big))
  for (j in 1:ncol(x.big)) {
    mylm <- lm(y ~ x.big[, j],
               weights = c(rep(N2, N1), rep(N1, N2)))
    res[1:2, j] <- mylm$coefficients
    res[3, j] <- summary(mylm)$r.squared
  }
  res
}

get_res2_class <- function(x.big, y) {
  rbind(CoeffsClass(x.big, y), RsqClass(x.big, y))
}

# In a case of classification
y <- c(rep(1, N1), rep(-1, N2))
y2 <- x + rnorm(length(x), 0, 0.1)
x.big <- as.big.matrix(cbind(x, x+1, 2*x))
x.big2 <- as.big.matrix(round(cbind(x, x+1, 2*x)))
x.big3 <- as.big.matrix(round(cbind(x, x+1, 2*x)), type = "char")

test_that("equality with lm in case of classification", {
  expect_equal(get_res2_class(x.big, y),
               get_res_weights(x.big, y))
  expect_equal(get_res2_class(x.big2, y),
               get_res_weights(x.big2, y))
  expect_equal(get_res2_class(x.big3, y),
               get_res_weights(x.big3, y))
})
test_that("expect error from regression", {
  expect_error(CoeffsReg(x.big, y))
  expect_error(RsqReg(x.big, y))
})

################################################################################

get_res <- function(x.big, y) {
  res <- matrix(0, 3, ncol(x.big))
  for (j in 1:ncol(x.big)) {
    mylm <- lm(y ~ x.big[, j])
    res[1:2, j] <- mylm$coefficients
    res[3, j] <- summary(mylm)$r.squared
  }
  res
}

get_res2_reg <- function(x.big, y) {
  rbind(CoeffsReg(x.big, y), RsqReg(x.big, y))
}

# In a case of regression
test_that("equality with lm in case of regression", {
  expect_equal(get_res2_reg(x.big, y2),
               get_res(x.big, y2))
  expect_equal(get_res2_reg(x.big2, y2),
               get_res(x.big2, y2))
  expect_equal(get_res2_reg(x.big3, y2),
               get_res(x.big3, y2))
})
test_that("expect error from classification", {
  expect_error(CoeffsClass(x.big, y2))
  expect_error(RsqClass(x.big, y2))
})

################################################################################

ind.train <- sort(sample(length(x), length(x) / 2))

get_res_train <- function(x.big, y) {
  res <- matrix(0, 3, ncol(x.big))
  for (j in 1:ncol(x.big)) {
    mylm <- lm(y[ind.train] ~ x.big[ind.train, j])
    res[1:2, j] <- mylm$coefficients
    res[3, j] <- summary(mylm)$r.squared
  }
  res
}

get_res2_reg_train <- function(x.big, y) {
  rbind(CoeffsReg(x.big, y, ind.train), RsqReg(x.big, y, ind.train))
}

# With only half of the data
test_that("equality with lm in case of regression with half of the data", {
  expect_equal(get_res2_reg_train(x.big, y2),
               get_res_train(x.big, y2))
  expect_equal(get_res2_reg_train(x.big2, y2),
               get_res_train(x.big2, y2))
  expect_equal(get_res2_reg_train(x.big3, y2),
               get_res_train(x.big3, y2))
})
test_that("expect error from classification", {
  expect_error(CoeffsClass(x.big, y2))
  expect_error(RsqClass(x.big, y2))
})

################################################################################

options(opt.save)

################################################################################
