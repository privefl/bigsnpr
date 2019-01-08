################################################################################

context("READ_BED")

test <- snp_attachExtdata()
G <- test$genotypes

bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
test2 <- gaston::read.bed.matrix(bedfile, verbose = FALSE)

################################################################################

test_that("same genotype matrix as gaston (reversed)", {
  expect_equivalent(G[], 2L - gaston::as.matrix(test2))
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

test_that("same sign as PLINK (no switch 0 <-> 2)", {
  plink <- download_plink()
  prefix <- sub("\\.bed$", "", bedfile)
  tmp <- tempfile()
  system(glue::glue("{plink} --bfile {prefix} --assoc --allow-no-sex --out {tmp}"))

  gwas <- big_univLogReg(G, test$fam$affection - 1L)
  sumstats <- bigreadr::fread2(paste0(tmp, ".assoc"))
  expect_gt(cor(gwas$estim, log(sumstats$OR)), 0.99)
})

################################################################################
