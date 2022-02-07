################################################################################

context("ATTACH")

################################################################################

test <- snp_attachExtdata()
expect_null(test$savedIn)
bkfile <- test$genotypes$backingfile
rdsfile <- sub_bk(bkfile, ".rds")

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

rds2bk <- function(rdsfile) sub("\\.rds$", ".bk", rdsfile)

rdsfile0 <- test_path("testdata", "before_readonly.rds")
rdsfile <- tempfile(fileext = ".rds")
file.copy(from = rdsfile0, to = rdsfile)
file.copy(from = rds2bk(rdsfile0), to = rds2bk(rdsfile))

test <- readRDS(rdsfile)
expect_false(exists("is_read_only", test$genotypes))

expect_message(test2 <- snp_attach(rdsfile),
               "[FBM from an old version? Reconstructing..|You should use `snp_save()`.]")
G2 <- test2$genotypes
expect_true(exists("is_read_only", G2))
G2[1] <- 3
expect_identical(G2[1], NA_real_)
G2$is_read_only <- TRUE
expect_error(G2[1] <- 4, "This FBM is read-only.")

################################################################################
