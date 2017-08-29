################################################################################

context("READ_BED")

test <- snp_attachExtdata()
G <- test$genotypes

bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
test2 <- gaston::read.bed.matrix(bedfile, verbose = FALSE)

################################################################################

test_that("same genotype matrix as gaston", {
  expect_equivalent(G[], gaston::as.matrix(test2))
  expect_equivalent(test$fam, test2@ped[1:6])
  expect_equivalent(test$map, test2@snps[1:6])
})

test_that("good class", {
  expect_s3_class(test, "bigSNP")
})

################################################################################

path <- sub("\\.bk$", "", G$backingfile)

test_that("Error: already exists", {
  expect_error(snp_readBed(bedfile, backingfile = path),
               sprintf("File '%s' already exists.", paste0(path, ".bk")),
               fixed = TRUE)
})

################################################################################
