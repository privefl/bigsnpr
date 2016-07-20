################################################################################

library(bigsnpr)
library(testthat)
context("READPLINK")

if (!dir.exists("backingfiles")) dir.create("backingfiles")
if (file.exists("backingfiles/test.bk")) file.remove("backingfiles/test.bk")

bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
test <- BedToBig(bedfile, 50, "test", "backingfiles")

# # requireNamespace("snpStats", quietly = TRUE)
# test2 <- snpStats::read.plink(bedfile)
# mat.test2 <- as(test2$genotypes, "numeric")
#
# ################################################################################
#
# test_that("same genotype matrix as snpStats", {
#   expect_equal(sum(abs(test$genotypes[,] - mat.test2)), 0)
# })
#
# test_that("good class", {
#   expect_match(class(test), "bigSNP")
# })
#
# test_that("data.frames of same size", {
#   expect_equal(dim(test$fam), dim(test2$fam))
#   expect_equal(dim(test$map), dim(test2$map))
# })

################################################################################

if (!dir.exists("backingfiles")) dir.create("backingfiles")
if (file.exists("backingfiles/test2.bk")) file.remove("backingfiles/test2.bk")

pedfile <- system.file("extdata", "example.ped.gz", package = "bigsnpr")

test3 <- PedToBig(pedfile, 50, "test2", "backingfiles")

test_that("Same types in data.frames", {
  expect_equal(sapply(test$fam, typeof), sapply(test3$fam, typeof))
  expect_equal(sapply(test$map, typeof), sapply(test3$map, typeof))
})

test_that("good class", {
  expect_match(class(test3), "bigSNP")
})

test_ref <- function(code) {
  check1 <- all(code %in% c(2, 20, 11, NA))
  check2 <- all(code %in% c(0, 11, 22, NA))

  check1 | check2
}

code.gen <- 10 * test$genotypes[,] + test3$genotypes[,]
check <- apply(code.gen, 2, test_ref)

test_that("Same genotypes depending on NA/ref", {
  expect_true(all(check))
})
