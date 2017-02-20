################################################################################

context("WRITEPLINK")

N <- round(runif(1, 10, 100))
M <- round(runif(1, 100, 1000))

fake <- snp_fake(N, M)
fake$genotypes[] <- sample(c(0:2, NA), size = N * M, replace = TRUE)

test_that("Error: already exists", {
  expect_error(BedToBig(bedfile, 50, "test_doc"))
})


################################################################################

# write the object as a bed/bim/fam object
tmpfile <- tempfile()
bed <- snp_writeBed(fake, paste0(tmpfile, ".bed"))

test_that("Error: already exists", {
  expect_error(snp_writeBed(fake, bed),
               sprintf("File %s already exists", bed), fixed = TRUE)
})


################################################################################

# read this new file for the first time
fake2 <- snp_attach(snp_readBed(bed, backingfile = basename(tmpfile),
                                backingpath = dirname(tmpfile)))

test_that("same content as written bigSNP", {
  expect_equal(fake$genotypes[,], fake2$genotypes[,])
  expect_equal(fake$fam, fake2$fam)
  expect_equal(fake$map, fake2$map)
})

################################################################################
