################################################################################

context("MODIF-SAVE")

if (!dir.exists("backingfiles")) dir.create("backingfiles")
if (file.exists("backingfiles/test4.bk")) file.remove("backingfiles/test4.bk")

bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
test <- BedToBig(bedfile, 50, "test4", "backingfiles")

################################################################################

# Just after reading
test_that("slot pop is initially NULL", {
  expect_null(test$fam$pop)
})

# Get populations from slot "family.ID"
test <- GetPops(test)
test_that("slot pop is now equal to slot family.ID", {
  expect_equal(test$fam$pop, test$fam$family.ID)
})

# Get populations clusters from external files
files <- system.file("extdata", paste0("cluster", 1:3), package = "bigsnpr")
test <- GetPops(x = test, pop.files = files,
                col.sample.ID = 2, col.family.ID = 3)
test_that("slot pop is now equal to slot 3 different pops", {
  expect_equal(rle(test$fam$pop)$values,
               c("POP-1-5", "POP-6-11", "POP-12-19"))
})

################################################################################
