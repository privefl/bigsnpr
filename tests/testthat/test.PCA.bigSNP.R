################################################################################

context("PCA")

if (!dir.exists("backingfiles")) dir.create("backingfiles")
if (file.exists("backingfiles/test5.bk")) file.remove("backingfiles/test5.bk")

bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
test <- BedToBig(bedfile, 50, "test5", "backingfiles")

################################################################################


# Get the first 10 PCs
test2 <- PCA.bigSNP(test, block.size = 1000, k = 10)

test_that("Sd of PCs is decreasing", {
  expect_lt(max(diff(apply(test2, 2, sd))), 0)
})

################################################################################

# Train PCs with only half the data
n <- nrow(test$genotypes)
ind.train <- sort(sample(n, n/2))
test3 <- PCA.bigSNP(test, block.size = 1000, k = 10,
                    ind.train = ind.train)

cor.2PC <- abs(diag(cor(test2, test3))[1:2])
test_that("Results correlated even with half the data", {
  expect_gt(min(cor.2PC), 0.9)
})

################################################################################

# Compare with prcomp
test4 <- prcomp(test$genotypes[,], scale. = TRUE)

cor.2PC2 <- abs(diag(cor(test2[, 1:2], test4$x[, 1:2])))
test_that("Results correlated with prcomp's result", {
  expect_gt(min(cor.2PC2), 0.9)
})

################################################################################
