################################################################################

context("BED_RANDOM_SVD")

options(bigstatsr.check.parallel.blas = FALSE)

skip_if(is_cran)

################################################################################

# No missing value -> with {bigstatsr}
bigSNP <- snp_attachExtdata()
G <- bigSNP$genotypes

test_that("Scales work", {
  sc1 <- snp_scaleBinom()(G)
  sc2 <- snp_scaleAlpha()(G)
  expect_equal(sc1, sc2)
  expect_equal(snp_scaleAlpha(alpha = 0)(G)$scale, rep(1, ncol(G)))
})

POS <- bigSNP$map$physical.pos
CHR0 <- sort(rep_len(1:2, ncol(G)))
ind.keep0 <- snp_clumping(G, infos.chr = CHR0, infos.pos = POS,
                          exclude = which(CHR0 == 1))
expect_true(all(CHR0[ind.keep0] == 2))

CHR <- bigSNP$map$chromosome
ind.keep <- snp_clumping(G, infos.chr = CHR, infos.pos = POS)
obj.svd <- big_randomSVD(G, fun.scaling = snp_scaleBinom(), ind.col = ind.keep)

# Try on bed file directly (still no missing value)
obj.bed <- bed(system.file("extdata", "example.bed", package = "bigsnpr"))
expect_error(bed_clumping(G))
ind.keep2 <- bed_clumping(obj.bed)
expect_identical(ind.keep2, ind.keep)
expect_gt(min(bed_clumping(obj.bed, exclude = 1:100)), 100)
expect_error(big_randomSVD(obj.bed, ind.col = ind.keep),
             "'X' is not of class 'FBM'.")
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

expect_identical(bed_clumping(obj.bed2, S = bed_MAF(obj.bed2)$mac), ind.keep3)

# Compute whole eigen decomposition
K <- bed_tcrossprodSelf(obj.bed2, ind.col = ind.keep2)
eig <- eigen(K[], symmetric = TRUE)
expect_equal(sqrt(eig$values[1:10]), obj.svd3$d)
expect_gt(mean(sqrt(colSums(cor(eig$vectors[, 1:10], obj.svd3$u)^2))), 0.999)

################################################################################

expect_identical(bed_clumping(obj.bed2, ncores = NCORES), ind.keep3)
expect_equal(bed_randomSVD(obj.bed2, ind.col = ind.keep2, ncores = NCORES),
             obj.svd3)

################################################################################

# Test class bed
expect_identical(capture.output(obj.bed3 <- print(obj.bed)),
                 "A 'bed' object with 517 samples and 4542 variants.")
expect_identical(obj.bed3, obj.bed)
expect_equal(length(obj.bed), prod(dim(obj.bed)))

# Test counts and scaling
ind.row <- sample(nrow(obj.bed), 300)
ind.col <- sample(ncol(obj.bed), 4000)
expect_identical(bed_counts(obj.bed2, ind.row, ind.col),
                 big_counts(G, ind.row, ind.col))
expect_identical(bed_counts(obj.bed2, ind.row, ind.col, ncores = 2),
                 big_counts(G, ind.row, ind.col))
expect_identical(bed_counts(obj.bed2, ind.row, ind.col, byrow = TRUE),
                 big_counts(G, ind.row, ind.col, byrow = TRUE))
expect_identical(bed_counts(obj.bed2, ind.row, ind.col, byrow = TRUE, ncores = 2),
                 big_counts(G, ind.row, ind.col, byrow = TRUE))

G_nona <- snp_attachExtdata()$genotypes
expect_identical(bed_counts(obj.bed, ind.row, ind.col),
                 big_counts(G_nona, ind.row, ind.col))
expect_identical(bed_counts(obj.bed, ind.row, ind.col, ncores = 2),
                 big_counts(G_nona, ind.row, ind.col))
expect_identical(bed_MAF(obj.bed, ind.row, ind.col)$N,
                 rep(length(ind.row), length(ind.col)))
expect_identical(bed_MAF(obj.bed, ind.row, ind.col)$maf,
                 snp_MAF(G_nona, ind.row, ind.col))

expect_identical(bed_MAF(obj.bed, ind.row, ind.col)$af,
                 bed_MAF(obj.bed, ind.row, ind.col)$ac / (2 * length(ind.row)))
expect_identical(bed_MAF(obj.bed2, ind.row, ind.col)$af,
                 bed_MAF(obj.bed2, ind.row, ind.col)$ac /
                   (2 * bed_MAF(obj.bed2, ind.row, ind.col)$N))
expect_true(all(bed_MAF(obj.bed2, ind.row, ind.col)$af <= 1))
expect_true(all(bed_MAF(obj.bed2, ind.row, ind.col)$ac <=
                  (2 * bed_MAF(obj.bed2, ind.row, ind.col)$N)))

ind1 <- sample(nrow(obj.bed2), 300)
ind2 <- sample(ncol(obj.bed2), 3000)
expect_identical(bed_scaleBinom(obj.bed2, ind1, ind2)$center,
                 2 * bed_MAF(obj.bed2, ind1, ind2, ncores = 2)$af)

################################################################################
