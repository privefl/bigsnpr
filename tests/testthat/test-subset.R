################################################################################

context("SUBSET")

test <- snp_attachExtdata()
# backing file
file <- test$genotypes@description$filename
dir  <- test$genotypes@description$dirname
# dimensions
n <- nrow(test$genotypes)
m <- ncol(test$genotypes)

################################################################################

for (k in 1:10) {
  # taking only some individuals and some SNPs at random
  ind.row <- sort(sample(n, size = sample(50, 1)))
  ind.col <- sort(sample(m, size = sample(200, 1)))
  test2 <- snp_attach(subset(test, ind.row, ind.col))

  # classes
  test_that("good classes", {
    expect_s3_class(test2, "bigSNP")
    expect_s4_class(test2$genotypes, "BM.code.descriptor")
    expect_s3_class(test2$fam, "data.frame")
    expect_s3_class(test2$map, "data.frame")
  })

  # dimensions
  test_that("Check new values", {
    expect_equal(attach.BM(test2$genotypes)[,],
                 attach.BM(test$genotypes)[ind.row, ind.col])
    expect_equivalent(test2$fam, test$fam[ind.row, ])
    expect_equivalent(test2$map, test$map[ind.col, ])
  })

  # new backing files
  test_that("backing files' names", {
    expect_match(test2$genotypes@description$dirname, dir)
    expect_equal(test2$genotypes@description$filename,
                 sub("\\.bk$", paste0("_sub", k, ".bk"), file))
  })
}

################################################################################

# removing the 100th and 120th individuals
ind.row.del <- c(100, 120)
test3 <- snp_attach(subset(test, -ind.row.del))

# classes
test_that("good classes", {
  expect_s3_class(test3, "bigSNP")
  expect_s4_class(test3$genotypes, "BM.code.descriptor")
  expect_s3_class(test3$fam, "data.frame")
  expect_s3_class(test3$map, "data.frame")
})

# dimensions
test_that("Check new values", {
  expect_equal(attach.BM(test3$genotypes)[,],
               attach.BM(test$genotypes)[-ind.row.del, ])
  expect_equivalent(test3$fam, test$fam[-ind.row.del, ])
  expect_equivalent(test3$map, test$map)
})

# new backing files
test_that("backing files' names", {
  expect_match(test3$genotypes@description$dirname, dir)
  expect_equal(test3$genotypes@description$filename,
               sub("\\.bk$", paste0("_sub", k + 1, ".bk"), file))
})

################################################################################
