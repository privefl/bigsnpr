################################################################################

context("FAKE")

################################################################################

n <- round(runif(1, 10, 100))
m <- round(runif(1, 10, 100))
test <- snp_fake(n, m)

test_that("good dimensions", {
  expect_equal(dim(test$genotypes), c(n, m))
  expect_equal(dim(test$fam), c(n, 6))
  expect_equal(dim(test$map), c(m, 6))
})

test_that("good classes", {
  expect_s3_class(test, "bigSNP")
  expect_s4_class(test$genotypes, "FBM.code256")
  expect_s3_class(test$fam, "data.frame")
  expect_s3_class(test$map, "data.frame")
})

################################################################################

G <- test$genotypes

test_that("only missing values", {
  expect_true(all(is.na(G[])))
})

# Modify the genotype `big.matrix`
G[] <- sample(as.raw(0:3), size = length(G), replace = TRUE)

counts <- big_counts(G)

test_that("counts 0, 1, 2, NA", {
  expect_equal(dim(counts), c(4, m))
})

################################################################################
