################################################################################

context("READPLINK")

if (!dir.exists("backingfiles")) dir.create("backingfiles")
if (file.exists("backingfiles/test.bk")) file.remove("backingfiles/test.bk")

bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
test <- BedToBig(bedfile, 50, "test", "backingfiles")

test2 <- gaston::read.bed.matrix(bedfile)

################################################################################

test_that("same genotype matrix as gaston", {
  expect_equivalent(test$genotypes[,], as.matrix(test2))
})

test_that("Error", {
  expect_error(BedToBig(bedfile, 50, "test", "back"))
  expect_error(BedToBig(bedfile, 50, "test", "backingfiles"))
})

test_that("good class", {
  expect_match(class(test), "bigSNP")
})

test_that("data.frames of same size as matrix", {
  expect_equal(dim(test$fam), c(nrow(test$genotypes), 6))
  expect_equal(dim(test$map), c(ncol(test$genotypes), 6))
})

################################################################################

# if (!dir.exists("backingfiles")) dir.create("backingfiles")
# if (file.exists("backingfiles/test2.bk")) file.remove("backingfiles/test2.bk")
#
# pedfile <- system.file("extdata", "example.ped.gz", package = "bigsnpr")
#
# test3 <- PedToBig(pedfile, 50, "test2", "backingfiles")
#
# test_that("Same types in data.frames", {
#   expect_equal(sapply(test$fam, typeof), sapply(test3$fam, typeof))
#   expect_equal(sapply(test$map, typeof), sapply(test3$map, typeof))
# })
#
# test_that("good class", {
#   expect_match(class(test3), "bigSNP")
# })
#
# test_ref <- function(code) {
#   check1 <- all(code %in% c(2, 20, 11, NA))
#   check2 <- all(code %in% c(0, 11, 22, NA))
#
#   check1 | check2
# }
#
# code.gen <- 10 * test$genotypes[,] + test3$genotypes[,]
# check <- apply(code.gen, 2, test_ref)
#
# test_that("Same genotypes depending on NA/ref", {
#   expect_true(all(check))
# })
