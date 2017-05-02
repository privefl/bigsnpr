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

# change the code
fake$genotypes@code[4] <- 2
# write the object as a bed/bim/fam object
bed <- snp_writeBed(fake, bedfile = tempfile(fileext = ".bed"))
tmpfile <- tempfile()
fake3 <- snp_attach(snp_readBed(bed, backingfile = basename(tmpfile),
                                backingpath = dirname(tmpfile)))

test_that("same content as written bigSNP (no more NAs)", {
  expect_equal(fake3$genotypes@code, CODE_012)
  expect_equal(sum(is.na(attach.BM(fake3$genotypes)[,])), 0)
  expect_equal(attach.BM(fake$genotypes)[,],
               attach.BM(fake3$genotypes)[,])
  expect_equal(fake$fam, fake3$fam)
  expect_equal(fake$map, fake3$map)
})

################################################################################

ind.row <- sample(N, 10)
ind.col <- sample(M, 50)
# write the object as a smaller bed/bim/fam object
bed <- snp_writeBed(fake, bedfile = tempfile(fileext = ".bed"),
                    ind.row = ind.row, ind.col = ind.col)
tmpfile <- tempfile()
fake4 <- snp_attach(snp_readBed(bed, backingfile = basename(tmpfile),
                                backingpath = dirname(tmpfile)))

test_that("same content as written bigSNP (with subsetting)", {
  expect_equal(fake4$genotypes@code, CODE_012)
  expect_equal(sum(is.na(attach.BM(fake3$genotypes)[,])), 0)
  expect_equal(attach.BM(fake$genotypes)[ind.row, ind.col],
               attach.BM(fake4$genotypes)[,])
  expect_equivalent(fake$fam[ind.row, ], fake4$fam)
  expect_equivalent(fake$map[ind.col, ], fake4$map)
})

################################################################################

rm(X)

################################################################################
