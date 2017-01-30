################################################################################

context("MODIF-SAVE")

test <- snp_attachExample("test-modif-save")

################################################################################

# Just after reading
test_that("slot pop is NULL", {
  expect_null(test$fam$pop)
})

# Get populations clusters from external files
files <- system.file("extdata", paste0("cluster", 1:3), package = "bigsnpr")
test <- snp_getPops(x = test, pop.files = files,
                    col.sample.ID = 2, col.family.ID = 3)
test_that("slot pop is now equal to slot 3 different pops", {
  expect_equal(rle(test$fam$family.ID)$values,
               c("POP-1-5", "POP-6-11", "POP-12-19"))
})

test_that("Warning if not all individuals are matched", {
  expect_warning(test <- snp_getPops(x = test, pop.files = files[1],
                                     col.sample.ID = 2, col.family.ID = 3))
})

################################################################################

# test$fam$affection <- sample(c(1, 2), size = nrow(test$fam), replace = TRUE)
#
# test_that("slot pheno is initially NULL", {
#   expect_null(test$fam$pheno)
# })
#
# test <- GetPhenos(test)
#
# test_that("slot pheno is in {-1, 1}", {
#   expect_equal(sort(unique(test$fam$pheno)), c(-1, 1))
# })
#
# test_that("expect warning from getting NAs", {
#   expect_warning(test <- GetPhenos(test, coded01 = TRUE))
# })
#
# test$fam$affection <- test$fam$affection - 1
# test <- GetPhenos(test, coded01 = TRUE)
#
# test_that("slot pheno is in {-1, 1}", {
#   expect_equal(sort(unique(test$fam$pheno)), c(-1, 1))
# })

################################################################################
