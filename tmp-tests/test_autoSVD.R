################################################################################

context("AUTO_SVD")

################################################################################

test <- snp_attachExtdata()
G   <- test$genotypes
CHR <- test$map$chromosome
POS <- test$map$physical.pos / 10
POS2 <- round(POS + 1)

################################################################################

expect_error(snp_autoSVD(G, CHR, POS), "only integers")

expect_is(obj.svd <- snp_autoSVD(G, CHR2), "big_SVD")
expect_identical(attr(obj.svd, "lrldr"), LD.wiki34[0, 1:3])

obj.svd2 <- snp_autoSVD(G, CHR, POS2, size = 5)
expect_gt(length(attr(obj.svd2, "subset")), length(attr(obj.svd, "subset")))
obj.svd3 <- snp_autoSVD(G, CHR, POS2, size = 5, is.size.in.bp = TRUE)
expect_lt(length(attr(obj.svd3, "subset")), length(attr(obj.svd2, "subset")))

obj.svd4 <- snp_autoSVD(G, CHR, roll.size = 0)
expect_lt(length(attr(obj.svd4, "subset")), length(attr(obj.svd, "subset")))

obj.svd5 <- snp_autoSVD(G, CHR, thr.r2 = 1, roll.size = 0)
expect_gt(min(diag(cor(obj.svd5$u, obj.svd$u))[1:3]^2), 0.98)

################################################################################
