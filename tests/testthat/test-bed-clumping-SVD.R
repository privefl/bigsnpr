################################################################################

# library(bigsnpr)
# library(testthat)
context("BED_RANDOM_SVD")

################################################################################

# No missing value -> with {bigstatsr}
bigSNP <- snp_attachExtdata()
G <- bigSNP$genotypes
CHR <- bigSNP$map$chromosome
POS <- bigSNP$map$physical.pos
ind.keep <- snp_clumping(G, infos.chr = CHR, infos.pos = POS)
obj.svd <- big_randomSVD(G, fun.scaling = snp_scaleBinom(), ind.col = ind.keep)

# Try on bed file directly (still no missing value)
obj.bed <- bed(system.file("extdata", "example.bed", package = "bigsnpr"))
expect_error(bed_clumping(G))
ind.keep2 <- bed_clumping(obj.bed)
expect_identical(ind.keep2, ind.keep)
expect_error(bed_randomSVD(G, ind.col = ind.keep),
             "'obj.bed' is not of class 'bed'.")
obj.svd2 <- bed_randomSVD(obj.bed, ind.col = ind.keep)
expect_equal(colMeans(obj.svd2$u), rep(0, 10))
expect_equal(obj.svd2, obj.svd)

# Try with missing value (supported for bed files only)
G[sample(length(G), length(G) / 20)] <- 3
obj.bed2 <- bed(snp_writeBed(bigSNP, tempfile(fileext = ".bed")))
ind.keep3 <- bed_clumping(obj.bed2)
expect_gt(length(intersect(ind.keep3, ind.keep2)) /
            length(union(ind.keep3, ind.keep2)), 0.95)
obj.svd3 <- bed_randomSVD(obj.bed2, ind.col = ind.keep2)
expect_gt(mean(sqrt(colSums(cor(obj.svd3$u, obj.svd2$u)^2))), 0.9)
expect_equal(obj.svd3$d, obj.svd2$d, tolerance = 0.1)
expect_true(all(obj.svd3$d < obj.svd2$d))

# Compute whole eigen decomposition
K <- bed_tcrossprodSelf(obj.bed2, ind.col = ind.keep2)
eig <- eigen(K[], symmetric = TRUE)
expect_equal(sqrt(eig$values[1:10]), obj.svd3$d)
expect_gt(mean(sqrt(colSums(cor(eig$vectors[, 1:10], obj.svd3$u)^2))), 0.999)

################################################################################

# https://github.com/hadley/testthat/issues/567
Sys.unsetenv("R_TESTS")

not_cran <- identical(Sys.getenv("BIGSNPR_CRAN"), "false")
NCORES <- `if`(not_cran, 2, 1)

expect_identical(bed_clumping(obj.bed2, ncores = NCORES), ind.keep3)
expect_equal(bed_randomSVD(obj.bed2, ind.col = ind.keep2, ncores = NCORES),
             obj.svd3)

################################################################################
