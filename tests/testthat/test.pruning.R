################################################################################

context("PRUNE")

if (!dir.exists("backingfiles")) dir.create("backingfiles")
if (file.exists("backingfiles/test3.bk")) file.remove("backingfiles/test3.bk")

bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
test <- BedToBig(bedfile, 50, "test3", "backingfiles")

test_that("Phenotypes needs to be defined", {
  expect_error(Prune(test))
})

################################################################################

test$fam$pheno <- sample(c(-1, 1), size = nrow(test$fam),
                         replace = TRUE)
ind <- Prune(test)
# ind2 <- Prune(test, ncores = 2)
ind3 <- Prune(test, thr.pvalue = 2)
ind4 <- Prune(test, thr.corr = 0.5)

################################################################################

# test_that("Same results with parallelism", {
#   expect_equal(ind, ind2)
# })

test_that("Lengths differ with thresholds", {
  expect_gte(length(ind), length(ind3))
  expect_lte(length(ind), length(ind4))
})

################################################################################
