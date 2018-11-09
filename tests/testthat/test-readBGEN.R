################################################################################

library(bigsnpr)
library(testthat)
context("READ_BGEN")

bgen_file <- system.file("testdata", "example.8bits.bgen", package = "bigsnpr")
variants <- readRDS(system.file("testdata", "bgen_variants.rds", package = "bigsnpr"))
dosages <- readRDS(system.file("testdata", "bgen_dosages.rds", package = "bigsnpr"))
IDs <- with(variants, paste(1, physical.pos, allele1, allele2, sep = "_"))
# variants 18 & 19 have identical IDs
excl <- c(18, 19)

################################################################################

test_that("raises some errors", {
  expect_error(snp_attach(snp_readBGEN(bgen_file, tempfile(), IDs)),
               "'list_snp_id' is not of class 'list'.", fixed = TRUE)
  expect_error(snp_attach(snp_readBGEN(bgen_file, tempfile(), list(c(IDs, "LOL")))),
               "Some variants have not been found")
})

test_that("same as package {rbgen}", {
  test <- snp_attach(snp_readBGEN(bgen_file, tempfile(), list(IDs)))
  G <- test$genotypes
  expect_identical(test$map[-excl, ], variants[-excl, ])
  expect_identical(G[, -excl][501], NA_real_)
  expect_equal(G[, -excl], round(dosages[, -excl], 2))
})

################################################################################

test_that("works with a subset of SNPs", {
  ind_snp <- setdiff(sample(length(IDs), 50), excl)
  test2 <- snp_attach(snp_readBGEN(bgen_file, tempfile(), list(IDs[ind_snp])))
  G2 <- test2$genotypes
  expect_equal(dim(G2), c(500, length(ind_snp)))
  expect_identical(test2$map, variants[ind_snp, ])
  expect_equal(G2[], round(dosages[, ind_snp], 2))
})

test_that("works with a subset of individuals", {
  ind_snp <- setdiff(sample(length(IDs), 50), excl)
  ind_row <- sample(500, 100)
  test3 <- snp_attach(
    snp_readBGEN(bgen_file, tempfile(), list(IDs[ind_snp]), ind_row))
  G3 <- test3$genotypes
  expect_equal(dim(G3), c(length(ind_row), length(ind_snp)))
  expect_identical(test3$map, variants[ind_snp, ])
  expect_equal(G3[], round(dosages[ind_row, ind_snp], 2))
})

################################################################################

test_that("works with multiple files", {
  ind_snp <- setdiff(sample(length(IDs), 50), excl)
  ind_row <- sample(500, 100)
  list_IDs <- split(IDs[ind_snp], sort(rep_len(1:3, length(ind_snp))))
  test4 <- snp_attach(
    snp_readBGEN(rep(bgen_file, 3), tempfile(), list_IDs, ind_row))
  G4 <- test4$genotypes
  expect_equal(dim(G4), c(length(ind_row), length(ind_snp)))
  expect_identical(test4$map, variants[ind_snp, ])
  expect_equal(G4[], round(dosages[ind_row, ind_snp], 2))
})

################################################################################
