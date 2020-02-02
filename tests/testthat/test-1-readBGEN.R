################################################################################

context("READ_BGEN")

skip_if_not_installed("dplyr")
skip_if_not_installed("dbplyr")
skip_if_not_installed("RSQLite")

################################################################################

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

ncores <- function() sample(1:2, 1)

################################################################################

test_that("raises some errors", {
  expect_error(snp_attach(snp_readBGEN(bgen_file, tempfile(), IDs, ncores = ncores())),
               "'list_snp_id' is not of class 'list'.", fixed = TRUE)
  expect_error(
    snp_attach(snp_readBGEN(bgen_file, tempfile(), list(c(IDs, "LOL")),
                            ncores = ncores())),
               "Wrong format of some variants.", fixed = TRUE)
})

################################################################################

test_that("format_snp_id() works as expected", {
  expect_error(format_snp_id(c("1_88169_C_T", "01_88169_C_T", "1:88169_C_T")),
               "Wrong format of some variants.", fixed = TRUE)
  expect_identical(format_snp_id(c("1_88169_C_T", "01_88169_C_T")),
                   c("01_88169_C_T", "01_88169_C_T"))
})

################################################################################

test_that("same as package {rbgen}", {
  test <- snp_attach(snp_readBGEN(bgen_file, tempfile(), list(IDs), ncores = ncores()))
  G <- test$genotypes
  expect_identical(test$map[-excl, ], variants[-excl, ])
  expect_identical(G[, -excl][501], NA_real_)
  expect_equal(G[, -excl], round(dosages[, -excl], 2))
})

################################################################################

setdiff2 <- function(x, y) x[!x %in% y]

test_that("works with a subset of SNPs", {
  ind_snp <- setdiff2(sample(length(IDs), 100, replace = TRUE), excl)
  expect_gte(length(ind_snp), 85)
  test2 <- snp_attach(snp_readBGEN(bgen_file, tempfile(), list(IDs[ind_snp]),
                                   ncores = ncores()))
  G2 <- test2$genotypes
  expect_equal(dim(G2), c(500, length(ind_snp)))
  expect_identical(test2$map, variants[ind_snp, ])
  expect_equal(G2[], round(dosages[, ind_snp], 2))
})

test_that("works with a subset of individuals", {
  ind_snp <- setdiff2(sample(length(IDs), 100, replace = TRUE), excl)
  ind_row <- sample(500, 100, replace = TRUE)
  test3 <- snp_attach(
    snp_readBGEN(bgen_file, tempfile(), list(IDs[ind_snp]), ind_row, ncores = ncores()))
  G3 <- test3$genotypes
  expect_equal(dim(G3), c(length(ind_row), length(ind_snp)))
  expect_identical(test3$map, variants[ind_snp, ])
  expect_equal(G3[], round(dosages[ind_row, ind_snp], 2))
})

################################################################################

test_that("works with multiple files", {
  ind_snp <- setdiff2(sample(length(IDs), 100, replace = TRUE), excl)
  ind_row <- sample(500, 100, replace = TRUE)
  list_IDs <- split(IDs[ind_snp], sort(rep_len(1:3, length(ind_snp))))
  test4 <- snp_attach(
    snp_readBGEN(rep(bgen_file, 3), tempfile(), list_IDs, ind_row, ncores = ncores()))
  G4 <- test4$genotypes
  expect_equal(dim(G4), c(length(ind_row), length(ind_snp)))
  expect_identical(test4$map, variants[ind_snp, ])
  expect_equal(G4[], round(dosages[ind_row, ind_snp], 2))
})

################################################################################

test_that("read as random hard calls", {
  G <- snp_attach(snp_readBGEN(bgen_file, tempfile(), list(IDs)))$genotypes[]
  all_G <- replicate(100, simplify = FALSE, {
    test <- snp_readBGEN(bgen_file, tempfile(), list(IDs), read_as = "random")
    snp_attach(test)$genotypes[]
  })
  expect_false(any(sapply(all_G[-1], identical, all_G[[1]])))
  G_mean <- Reduce('+', all_G) / length(all_G)
  # G[1:5, 1:5] - G_mean[1:5, 1:5]
  # plot(G_mean, G, pch = 20, col = scales::alpha("black", 0.1))
  expect_identical(G_mean[is.na(G)], rep(NA_real_, 2))
  expect_equal(G_mean[!is.na(G)], G[!is.na(G)], tolerance = 0.1)
  lm_coef <- stats::lm(G_mean[!is.na(G)] ~ G[!is.na(G)])$coef
  expect_equal(unname(lm_coef), c(0, 1), tolerance = 1e-3)
})

################################################################################
