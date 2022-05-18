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

  simu <- snp_simuPheno(G, 0.2, 500, alpha = -2)
  sd <- sqrt(big_colstats(G, ind.col = simu$set)$var)
  cor1 <- pcor(simu$allelic_effects^2, sd, z = NULL)
  expect_true(all(cor1 < 0))

  simu2 <- snp_simuPheno(G, 0.2, 500, alpha = 0)
  sd2 <- sqrt(big_colstats(G, ind.col = simu2$set)$var)
  cor2 <- pcor(simu2$allelic_effects^2, sd2, z = NULL, alpha = 1e-5)
  expect_true(cor2[2] < 0 && cor2[3] > 0)
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
