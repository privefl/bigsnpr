################################################################################

context("SIMU_PHENO")

################################################################################

test_that("snp_simuPheno() works", {

  skip_if(is_cran)

  bigsnp <- snp_attachExtdata()
  G <- bigsnp$genotypes

  ind.row <- sample(nrow(G), 50)
  simu <- snp_simuPheno(G, 0.2, 20, ind.row = ind.row, ncores = 2)
  expect_length(simu$pheno, 50)
  expect_length(simu$set, 20)
  expect_length(simu$effects, 20)
  var_g <- drop(var(scale(G[ind.row, simu$set]) %*% simu$effects))
  expect_equal(var_g, 0.2)

  ind.possible <- sample(ncol(G), 50)
  simu <- snp_simuPheno(G, 0.2, 20, ind.possible = ind.possible, ncores = 2)
  expect_length(simu$set, 20)
  expect_length(simu$effects, 20)
  expect_true(all(simu$set %in% ind.possible))

  all_cor <- replicate(20, {
    h2 <- runif(1, min = 0.2, max = 0.8)
    M <- sample(ncol(G), 1)
    simu <- snp_simuPheno(G, h2, M)
    var_g <- drop(var(scale(G[, simu$set]) %*% simu$effects))
    expect_equal(var_g, h2)

    gwas <- big_univLinReg(G, simu$pheno)
    cor(gwas$estim[simu$set], simu$effects)
  })
  # Not sure where there can be NAs, but it happened once
  expect_gt(median(all_cor, na.rm = TRUE), 0.15)
})

################################################################################

test_that("alpha in snp_simuPheno() works", {

  skip_if(is_cran)

  bigsnp <- snp_attachExtdata()
  G <- bigsnp$genotypes

  FUN <- function(x, log_var, beta2) {
    S <- 1 + x[[1]]; sigma2 <- x[[2]]
    S * sum(log_var) + length(log_var) * log(sigma2) + sum(beta2 / exp(S * log_var)) / sigma2
  }

  DER <- function(x, log_var, beta2) {
    S <- 1 + x[[1]]; sigma2 <- x[[2]]
    res1 <- sum(log_var) - sum(log_var * beta2 / exp(S * log_var)) / sigma2
    res2 <- length(log_var) / sigma2 - sum(beta2 / exp(S * log_var)) / sigma2^2
    c(res1, res2)
  }

  SEQ <- seq(-1.8, 0.8, by = 0.05)
  res <- sapply(SEQ, function(alpha) {
    simu <- snp_simuPheno(G, 0.2, 500, alpha = alpha)
    log_var <- log(big_colstats(G, ind.col = simu$set)$var)
    beta2 <- simu$effects^2
    optim(par = c(0.2 / 500, 0.5), fn = FUN, gr = DER, method = "L-BFGS-B",
          lower = c(-1.5, 0.2 / 5000), upper = c(0.5, 0.2 / 50),
          log_var = log_var, beta2 = beta2)$par[1]
  })
  # plot(res, SEQ); abline(0, 1, col = "red", lwd = 2)
  expect_equal(res, SEQ, tolerance = 0.5)
  expect_true(all(res <= 0.5))
  expect_true(all(res >= -1.5))
})

################################################################################

test_that("prob in snp_simuPheno() works", {

  bigsnp <- snp_attachExtdata()
  G <- bigsnp$genotypes

  simu <- snp_simuPheno(G, 0.2, 500,
                        prob = c(rep(10, 100), rep(0, 100), rep(1, ncol(G) - 200)))
  expect_gt(sum(simu$set %in% 1:100), 20)  ## 11 are expected if random
  expect_equal(sum(simu$set %in% 101:200), 0)
})

################################################################################
