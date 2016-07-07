################################################################################

skip_on_cran()
context("READPLINK")

if (!file.exists("backingfiles")) dir.create("backingfiles")
if (file.exists("backingfiles/test")) file.remove("backingfiles/test")
if (file.exists("backingfiles/test.desc")) file.remove("backingfiles/test.desc")
if (file.exists("backingfiles/test.rds")) file.remove("backingfiles/test.rds")

bedfile <- system.file("extdata", "example.bed", package = "mypack")
test <- BedToBig(bedfile, 50, "test", "backingfiles")

# requireNamespace("snpStats", quietly = TRUE)
test2 <- snpStats::read.plink(bedfile)
mat.test2 <- as(test2$genotypes, "numeric")

################################################################################

test_that("same genotype matrix as snpStats", {
  expect_equal(sum(abs(test$genotypes[,] - mat.test2)), 0)
})

test_that("data.frames of same size", {
  expect_equal(dim(test$fam), dim(test2$fam))
  expect_equal(dim(test$map), dim(test2$map))
})

################################################################################
