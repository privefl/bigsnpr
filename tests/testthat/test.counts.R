################################################################################

context("COUNTS")

# constructing a fake genotype big.matrix
a <- big.matrix(10, 8, type = "char", shared = FALSE)
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

res2 <- numeric(10)
res2[] <- c(2, 2, 1, 2, 1, 3, 2, 2, 3, 3)

################################################################################

# constructing a fake incomplete bigSNP with 10 individuals and 15 SNPs
# where the 5 first individuals are cases and the 5 last are controls.
fake <- list()
class(fake) <- "bigSNP"
fake$genotypes <- a
fake$fam$pheno <- c(rep(1, 5), rep(-1, 5))

################################################################################

# Get counts
test <- Counts(fake)
test_that("Same results as by hand", {
  expect_equal(max(abs(test$counts.col - res)), 0)
  expect_equal(test$counts.row, res2)
})

################################################################################
