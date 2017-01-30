################################################################################

context("READPLINK")

test <- snp_attachExample()

bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
test2 <- gaston::read.bed.matrix(bedfile)

################################################################################

test_that("same genotype matrix as gaston", {
  expect_equivalent(test$genotypes[,], as.matrix(test2))
})

test_that("Error: already exists", {
  expect_error(BedToBig(bedfile, 50, "test_doc"))
})

test_that("good class", {
  expect_match(class(test), "bigSNP")
})

test_that("data.frames of same size as matrix", {
  expect_equal(dim(test$fam), c(nrow(test$genotypes), 6))
  expect_equal(dim(test$map), c(ncol(test$genotypes), 6))
})

################################################################################
