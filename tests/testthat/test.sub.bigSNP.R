################################################################################

context("SUB")

if (!dir.exists("backingfiles")) dir.create("backingfiles")
if (file.exists("backingfiles/test6.bk")) file.remove("backingfiles/test6.bk")
if (file.exists("backingfiles/test6_sub1.bk"))
  file.remove("backingfiles/test6_sub1.bk")
if (file.exists("backingfiles/test6_sub1.desc"))
  file.remove("backingfiles/test6_sub1.desc")
if (file.exists("backingfiles/test6_sub2.bk"))
  file.remove("backingfiles/test6_sub2.bk")
if (file.exists("backingfiles/test6_sub2.desc"))
  file.remove("backingfiles/test6_sub2.desc")

bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
test <- BedToBig(bedfile, 50, "test6", "backingfiles")

################################################################################

# dimensions
n <- nrow(test$genotypes)
m <- ncol(test$genotypes)

# taking only the first 50 individuals and 500 SNPs at random
ind.row <- 1:50
ind.col <- sort(sample(ncol(test$genotypes), 500))
test2 <- sub.bigSNP(test, ind.row, ind.col)

# dimensions
test_that("Check new dimensions", {
  expect_equal(dim(test2$genotypes), c(50, 500))
  expect_equal(dim(test2$fam), c(50, 6))
  expect_equal(dim(test2$map), c(500, 6))
})

# new backing files
test_that("backing files' names", {
  expect_match(test2$backingpath, test$backingpath)
  expect_equal(test2$backingfile, paste0(test$backingfile, "_sub1"))
})

################################################################################

# removing the 100th and 120th individuals
ind.row.del <- c(100, 120)
test3 <- sub.bigSNP(test, -ind.row.del)

# dimensions
test_that("Check new dimensions", {
  expect_equal(dim(test3$genotypes), c(n-2, m))
  expect_equal(dim(test3$fam), c(n-2, 6))
  expect_equal(dim(test3$map), c(m, 6))
})

# new backing files
test_that("backing files' names", {
  expect_match(test3$backingpath, test$backingpath)
  expect_equal(test3$backingfile, paste0(test$backingfile, "_sub2"))
})

################################################################################
