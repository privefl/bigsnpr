################################################################################

context("PCA_PROJECT")

# test_that()

################################################################################

obj.bed <- bed(system.file("extdata", "example-missing.bed", package = "bigsnpr"))
ind.row <- sample(nrow(obj.bed), 100)
ind.col <- which(bed_MAF(obj.bed, ind.row)$mac > 5)
obj.svd <- bed_randomSVD(obj.bed, ind.row = ind.row, ind.col = ind.col)

ind.test <- setdiff(rows_along(obj.bed), ind.row)
expect_error(bed_projectSelfPCA(obj.svd, obj.bed, ind.row = ind.test),
             "'ind.col' can't be `NULL`.")
expect_error(bed_projectSelfPCA(obj.svd, obj.bed, ind.row = ind.test, ind.col = ind.col[-1]),
             "Incompatibility between dimensions.")
proj <- bed_projectSelfPCA(obj.svd, obj.bed,
                           ind.row = rows_along(obj.bed),
                           ind.col = ind.col)
expect_equal(proj$simple_proj[ind.row, ], predict(obj.svd), tolerance = 1e-4)

proj2 <- bed_projectPCA(obj.bed, obj.bed,
                        ind.row.new = ind.test,
                        ind.row.ref = ind.row,
                        strand_flip = FALSE,
                        roll.size = 5,
                        thr.r2 = 0.8,
                        verbose = FALSE)

obj.svd2 <- bed_autoSVD(obj.bed, ind.row = ind.row,
                        roll.size = 5, thr.r2 = 0.8, verbose = FALSE)
proj3 <- bed_projectSelfPCA(obj.svd2, obj.bed, ind.row = ind.test)
expect_equal(proj2, proj3)

################################################################################

obj.bed <- bed(system.file("extdata", "example.bed", package = "bigsnpr"))
ind.row <- sample(nrow(obj.bed), 400)
obj.svd <- bed_randomSVD(obj.bed, ind.row = ind.row)

ind.test <- setdiff(rows_along(obj.bed), ind.row)
expect_error(bed_projectSelfPCA(obj.svd, obj.bed, ind.row = ind.test),
             "'ind.col' can't be `NULL`.")
expect_error(bed_projectSelfPCA(obj.svd, obj.bed, ind.row = ind.test, ind.col = 1:5),
             "Incompatibility between dimensions.")
proj <- bed_projectSelfPCA(obj.svd, obj.bed,
                           ind.row = rows_along(obj.bed),
                           ind.col = cols_along(obj.bed))
expect_equal(proj$simple_proj[ind.row, ], predict(obj.svd), tolerance = 1e-4)

pop <- rep(1:3, c(143, 167, 207))
geom_med <- function(U) suppressWarnings(bigutilsr::geometric_median(U))
ref   <- unlist(by(predict(obj.svd)[, 2:3],         pop[ind.row],  geom_med))
pred1 <- unlist(by(proj$simple_proj[ind.test, 2:3], pop[ind.test], geom_med))
pred2 <- unlist(by(proj$OADP_proj[ind.test, 2:3],   pop[ind.test], geom_med))
expect_gt(sum(ref^2), sum(pred1^2))
expect_lt(sum((ref - pred2)^2), sum((ref - pred1)^2))

G <- snp_attach(snp_readBed(obj.bed$bedfile, backingfile = tempfile()))$genotypes
proj.2 <- snp_projectSelfPCA(obj.svd, G, ind.row = ind.test,
                             ind.col = cols_along(G), ncores = 2)
expect_identical(proj.2$obj.svd.ref, proj$obj.svd.ref)
expect_equal(proj.2$simple_proj, proj$simple_proj[ind.test, ])
expect_equal(proj.2$OADP_proj,   proj$OADP_proj[ind.test, ])


proj2 <- bed_projectPCA(obj.bed, obj.bed,
                        ind.row.new = ind.test,
                        ind.row.ref = ind.row,
                        strand_flip = FALSE,
                        roll.size = 10,
                        thr.r2 = 0.8,
                        verbose = FALSE)

obj.svd2 <- bed_autoSVD(obj.bed, ind.row = ind.row,
                        roll.size = 10, thr.r2 = 0.8, verbose = FALSE)
proj3 <- bed_projectSelfPCA(obj.svd2, obj.bed, ind.row = ind.test)
expect_equal(proj2, proj3, tolerance = 1e-6)

################################################################################
