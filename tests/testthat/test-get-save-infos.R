################################################################################

context("GET_SAVE_INFOS")

test <- snp_attachExtdata()

################################################################################

# Just after reading
test_that("family.ID with 19 pops after reading", {
  expect_equal(length(rle(test$fam$family.ID)$values), 19)
})

# Get populations clusters from external files
files <- system.file("extdata", paste0("cluster", 1:3), package = "bigsnpr")
infos <- snp_getSampleInfos(test, files, header = FALSE)[[1]]
test_that("slot pop is now equal to slot 3 different pops", {
  expect_equal(rle(infos)$values,
               c("POP-1-5", "POP-6-11", "POP-12-19"))
})

test_that("Warning if some individuals are not matched", {
  expect_warning(snp_getSampleInfos(test, files[1:2], header = FALSE))
})

################################################################################

# change slot
test$fam$family.ID <- infos
# re-attach
test <- snp_attach(sub("\\.bk$", ".rds", test$genotypes$backingfile))

test_that("family.ID with 19 pops after re-attaching", {
  expect_equal(length(rle(test$fam$family.ID)$values), 19)
})

################################################################################

# change slot
test$fam$family.ID <- infos
# save modifs
test <- snp_save(test)
# re-attach
test <- snp_attach(sub("\\.bk$", ".rds", test$genotypes$backingfile))

test_that("slot pop is now equal to slot 3 different pops", {
  expect_equal(rle(test$fam$family.ID)$values,
               c("POP-1-5", "POP-6-11", "POP-12-19"))
})

################################################################################
