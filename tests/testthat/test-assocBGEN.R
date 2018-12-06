################################################################################

context("ASSOC_BGEN")

# need to write bgen/bgi files because can't have binary files..
library(magrittr)
bgen_file <- tempfile(fileext = ".bgen")
system.file("testdata", "bgen_example.rds", package = "bigsnpr") %>%
  readRDS() %>% writeBin(bgen_file, useBytes = TRUE)
system.file("testdata", "bgi_example.rds",  package = "bigsnpr") %>%
  readRDS() %>% writeBin(paste0(bgen_file, ".bgi"), useBytes = TRUE)

variants <- readRDS(system.file("testdata", "bgen_variants.rds", package = "bigsnpr"))
dosages <- readRDS(system.file("testdata", "bgen_dosages.rds", package = "bigsnpr"))
IDs <- with(variants, paste(1, physical.pos, allele1, allele2, sep = "_"))
# variants 18 & 19 have identical IDs
excl <- c(18, 19)
excl2 <- union(which(is.na(dosages), arr.ind = TRUE)[, "col"], excl)

ncores <- function() sample(1:2, 1)
y <- rnorm(nrow(dosages))
ind <- rows_along(dosages)

################################################################################

test_that("raises some errors", {
  expect_error(snp_assocBGEN(bgen_file, IDs, y, ind, ncores = ncores()),
               "'list_snp_id' is not of class 'list'.", fixed = TRUE)
  expect_error(
    snp_assocBGEN(bgen_file, list(c(IDs, "LOL")), y, ind, ncores = ncores()),
    "Wrong format of some SNPs.", fixed = TRUE)
})

################################################################################

test_that("same as package {rbgen}", {
  test <- snp_assocBGEN(bgen_file, list(IDs), y, ind, ncores = ncores())
  lpval <- predict(big_univLinReg(as_FBM(dosages[, -excl2]), y))
  expect_equal(test[-excl2], lpval, tolerance = 1e-2)
})

################################################################################

test_that("works with a subset of SNPs", {
  ind_snp <- setdiff(sample(length(IDs), 50), excl2)
  test2 <- snp_assocBGEN(bgen_file, list(IDs[ind_snp]), y, ind, ncores = ncores())
  lpval2 <- predict(big_univLinReg(as_FBM(dosages[, ind_snp]), y))
  expect_equal(test2, lpval2, tolerance = 1e-2)
})

test_that("works with a subset of individuals", {
  ind_snp <- setdiff(sample(length(IDs), 50), excl2)
  ind_row <- sort(sample(500, 200))
  expect_error(
    snp_assocBGEN(bgen_file, list(IDs[ind_snp]), y, ind_row, ncores = ncores()),
    "Incompatibility between dimensions.", fixed = TRUE)
  test3 <- snp_assocBGEN(bgen_file, list(IDs[ind_snp]), y[ind_row], ind_row,
                         ncores = ncores())
  lpval3 <- predict(big_univLinReg(as_FBM(dosages[ind_row, ind_snp]), y[ind_row]))
  expect_equal(test3, lpval3, tolerance = 1e-2)
})

################################################################################

test_that("works with multiple files", {
  ind_snp <- setdiff(sample(length(IDs), 50), excl2)
  ind_row <- sort(sample(500, 200))
  list_IDs <- split(IDs[ind_snp], sort(rep_len(1:3, length(ind_snp))))
  expect_error(
    snp_assocBGEN(bgen_file, list_IDs, y[ind_row], ind_row, ncores = ncores()),
    "Incompatibility between dimensions.", fixed = TRUE)
  test4 <- snp_assocBGEN(rep(bgen_file, 3), list_IDs, y[ind_row], ind_row,
                         ncores = ncores())
  lpval4 <- predict(big_univLinReg(as_FBM(dosages[ind_row, ind_snp]), y[ind_row]))
  expect_equal(test4, lpval4, tolerance = 1e-2)
})

################################################################################
