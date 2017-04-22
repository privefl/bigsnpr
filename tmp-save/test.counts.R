################################################################################

context("COUNTS")

opt.save <- options(bigmemory.typecast.warning = FALSE,
                    bigmemory.default.shared = FALSE)

# constructing a fake genotype big.matrix
a <- big.matrix(10, 8, type = "char")
a[] <- c(NA, 1,  0,  1,  0,     2,  NA, NA, NA, NA,
         1,  1,  0,  2,  1,     2,  1,  NA, NA, 1,
         NA, NA, 0,  2,  NA,    NA, NA, 2,  1,  2,
         0,  1,  1,  2,  1,     NA, 2,  2,  1,  2,
         0,  1,  NA, NA, 0,     NA, 2,  2,  2,  NA,
         2,  2,  1,  1,  1,     0,  2,  0,  NA, 2,
         0,  NA, 2,  1,  1,     1,  0,  2,  2,  2,
         0,  2,  2,  NA, 0,     2,  1,  1,  0,  NA)

res <- matrix(0, 6, 8)
res[] <- c(2, 2, 0, 0, 0, 1,
           1, 3, 1, 0, 2, 1,
           1, 0, 1, 0, 1, 2,
           1, 3, 1, 0, 1, 3,
           2, 1, 0, 0, 0, 3,
           0, 3, 2, 2, 0, 2,
           1, 2, 1, 1, 1, 3,
           2, 0, 2, 1, 2, 1)

res2 <- c(2, 2, 1, 2, 1, 3, 2, 2, 3, 3)

################################################################################

# constructing a fake incomplete bigSNP with 10 individuals and 15 SNPs
# where the 5 first individuals are cases and the 5 last are controls.
fake <- list()
class(fake) <- "bigSNP"
fake$genotypes <- a
fake$fam$affection <- c(rep(1, 5), rep(-1, 5))

################################################################################

# Get counts
test <- snp_counts(fake)
test_that("Same results as by hand", {
  expect_equal(test$cols.cases, res[1:3, ])
  expect_equal(test$cols.controls, res[4:6, ])
  expect_equal(test$rows, res2)
})

################################################################################

# Get counts without phenos
test2 <- snp_counts(fake, has.pheno = FALSE)
test$cols <- test$cols.controls + test$cols.cases
test_that("Same results as by hand with pheno = FALSE", {
  expect_equal(test$cols, test2$cols)
  expect_equal(test2$rows, res2)
})

################################################################################

options(opt.save)

################################################################################
