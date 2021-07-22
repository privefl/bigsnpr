################################################################################

context("PROD_BGEN")

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
dosages  <- readRDS(system.file("testdata", "bgen_dosages.rds",  package = "bigsnpr"))
snp_info <- readRDS(system.file("testdata", "bgen_varinfo.rds",  package = "bigsnpr"))
IDs <- with(variants, paste(1, physical.pos, allele1, allele2, sep = "_"))
# variants 18 & 19 have identical IDs
excl <- c(18, 19)

ncores <- function() sample(1:2, 1)

################################################################################

test_that("same as with intermediate FBM", {

  expect_equal(rowSums(array(1:12, c(2, 2, 3)), dims = 2),
               matrix(c(15, 18, 21, 24), 2))

  G <- snp_readBGEN(bgen_file, tempfile(), list(IDs), ncores = ncores()) %>%
    snp_attach() %>%
    .$genotypes
  beta <- matrix(rnorm(ncol(G) * 3), ncol = 3)
  prod1 <- snp_prodBGEN(bgen_file, beta, list(IDs),
                        block_size = 2, ncores = ncores())
  expect_equal(prod1, G[] %*% beta, tolerance = 0.1)

  prod2 <- snp_prodBGEN(bgen_file, beta[, 1], list(IDs), ncores = ncores())
  expect_null(dim(prod2))
  expect_equal(prod2, prod1[, 1])

  prod3 <- snp_prodBGEN(bgen_file, beta[, 1, drop = FALSE], list(IDs),
                        ncores = ncores())
  expect_equal(ncol(prod3), 1)
  expect_equal(prod3, prod1[, 1, drop = FALSE])

  set.seed(1)
  prod3 <- snp_prodBGEN(bgen_file, beta, list(IDs), read_as = "random")
  set.seed(1)
  G2 <- snp_readBGEN(bgen_file, tempfile(), list(IDs), read_as = "random") %>%
    snp_attach() %>%
    .$genotypes
  expect_equal(prod3, G2[] %*% beta)
})

################################################################################

setdiff2 <- function(x, y) x[!x %in% y]

test_that("works with a subset of SNPs", {
  ind_snp <- setdiff2(sample(length(IDs), 50, replace = TRUE), excl)
  test2 <- snp_attach(snp_readBGEN(bgen_file, tempfile(), list(IDs[ind_snp]),
                                   ncores = ncores()))
  G2 <- test2$genotypes
  beta <- matrix(rnorm(length(ind_snp) * 3), ncol = 3)
  prod <- snp_prodBGEN(bgen_file, beta, list(IDs[ind_snp]),
                       block_size = 2, ncores = ncores())
  expect_equal(prod, G2[] %*% beta, tolerance = 0.1)
})

test_that("works with a subset of individuals", {
  ind_snp <- setdiff2(sample(length(IDs), 100, replace = TRUE), excl)
  ind_row <- sample(500, 100, replace = TRUE)
  test3 <- snp_attach(
    snp_readBGEN(bgen_file, tempfile(), list(IDs[ind_snp]), ind_row,
                 ncores = ncores()))
  G3 <- test3$genotypes
  beta <- matrix(rnorm(length(ind_snp) * 3), ncol = 3)
  prod <- snp_prodBGEN(bgen_file, beta, list(IDs[ind_snp]), ind_row,
                       block_size = 2, ncores = ncores())
  expect_equal(prod, G3[] %*% beta, tolerance = 0.1)
})

################################################################################

test_that("works with multiple files", {
  ind_snp <- setdiff2(sample(length(IDs), 50, replace = TRUE), excl)
  ind_row <- sample(500, 100, replace = TRUE)
  list_IDs <- split(IDs[ind_snp], sort(rep_len(1:3, length(ind_snp))))
  test4 <- snp_attach(
    snp_readBGEN(rep(bgen_file, 3), tempfile(), list_IDs, ind_row,
                 ncores = ncores()))
  G4 <- test4$genotypes
  beta <- matrix(rnorm(length(ind_snp) * 3), ncol = 3)
  prod <- snp_prodBGEN(rep(bgen_file, 3), beta, list_IDs, ind_row,
                       block_size = 2, ncores = ncores())
  expect_equal(prod, G4[] %*% beta, tolerance = 0.1)
})

################################################################################
