################################################################################

context("WRITE_BED")

N <- round(runif(1, 10, 100))
M <- round(runif(1, 100, 1000))

fake <- snp_fake(N, M)
X <- attach.BM(fake$genotypes)
X[] <- sample(as.raw(0:3), size = N * M, replace = TRUE)

test_that("Error: already exists", {
  expect_error(BedToBig(bedfile, 50, "test_doc"))
})

################################################################################

test_that("Write signed but read unsigned", {
  tmpfile <- tempfile()
  x <- as.raw(0:255)
  bigsnpr:::testWrite(x, tmpfile)
  test <- readBin(tmpfile, what = "raw", n = 256)
  expect_equal(test, x)
})

################################################################################

# write the object as a bed/bim/fam object
bed <- snp_writeBed(fake, bedfile = tempfile(fileext = ".bed"))

test_that("Error: already exists", {
  expect_error(snp_writeBed(fake, bedfile = bed),
               sprintf("File '%s' already exists", bed), fixed = TRUE)
})

################################################################################

# read this new file for the first time
tmpfile <- tempfile()
fake2 <- snp_attach(snp_readBed(bed, backingfile = basename(tmpfile),
                                backingpath = dirname(tmpfile)))

test_that("same content as written bigSNP", {
  expect_equal(attach.BM(fake$genotypes)[,],
               attach.BM(fake2$genotypes)[,])
  expect_equal(fake$fam, fake2$fam)
  expect_equal(fake$map, fake2$map)
})

################################################################################

rm(X)

################################################################################
