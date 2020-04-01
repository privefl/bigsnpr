################################################################################

context("SIMU_PHENO")

################################################################################

test_that("snp_simuPheno() works", {

  bigsnp <- snp_attachExtdata()
  G <- bigsnp$genotypes

  ind.row <- sample(nrow(G), 50)
  simu <- snp_simuPheno(G, 0.2, 20, ind.row = ind.row)
  expect_length(simu$pheno, 50)
  expect_length(simu$set, 20)
  expect_length(simu$effects, 20)
  var_g <- drop(var(scale(G[ind.row, simu$set]) %*% simu$effects))
  expect_equal(var_g, 0.2)

  ind.possible <- sample(ncol(G), 50)
  simu <- snp_simuPheno(G, 0.2, 20, ind.possible = ind.possible)
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
  expect_gt(median(all_cor), 0.15)
})

################################################################################
