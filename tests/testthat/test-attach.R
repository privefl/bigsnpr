################################################################################

context("ATTACH")

################################################################################

test <- snp_attachExtdata()
rdsfile <- test$savedIn
bkfile <- sub("\\.rds$", ".bk", rdsfile)

################################################################################

rdsfile.copy <- tempfile(fileext = ".rds")
expect_true(file.copy(rdsfile, rdsfile.copy))

test2 <- readRDS(rdsfile.copy)
expect_equal(test2$savedIn, test$savedIn)

bkfile.copy <- sub("\\.rds$", ".bk", rdsfile.copy)
expect_error(snp_attach(rdsfile.copy),
             sprintf("File '%s' doesn't exist", bkfile.copy))

expect_true(file.copy(bkfile, bkfile.copy))
test3 <- snp_attach(rdsfile.copy)
expect_equal(test3$savedIn, rdsfile.copy)
expect_equal(test3$genotypes@description$dirname, dirname(rdsfile.copy))
expect_s4_class(attach.BM(test3$genotypes), "BM.code")

################################################################################
