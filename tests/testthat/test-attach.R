################################################################################

context("ATTACH")

################################################################################

test <- snp_attachExtdata()
expect_null(test$savedIn)
bkfile <- test$genotypes$backingfile
rdsfile <- sub("\\.bk$", ".rds", bkfile)

################################################################################

rdsfile.copy <- tempfile(fileext = ".rds")
expect_true(file.copy(rdsfile, rdsfile.copy))

test2 <- readRDS(rdsfile.copy)
expect_null(test2$savedIn)

bkfile.copy <- sub("\\.rds$", ".bk", rdsfile.copy)
expect_error(snp_attach(rdsfile.copy))  # File doesn't exist

expect_true(file.copy(bkfile, bkfile.copy))
test3 <- snp_attach(rdsfile.copy)
expect_null(test3$savedIn)
expect_equal(test3$genotypes$backingfile, normalizePath(bkfile.copy))
expect_s4_class(test3$genotypes, "FBM.code256")

################################################################################
