################################################################################

context("READ_BED")

################################################################################

test_that("sub_bed() works", {
  expect_identical(sub_bed("toto.bed"), "toto")
  expect_identical(sub_bed("toto.bed", ".bim"), "toto.bim")
  expect_error(sub_bed("toto.bed2"),
               "Path 'toto.bed2' must have 'bed' extension.")
  expect_error(sub_bed("toto.bed", "bim"), "extension starting with '.'")
  expect_identical(sub_bed("toto.bed", "bim", stop_if_not_ext = FALSE), "totobim")
})

################################################################################

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

path <- sub_bk(G$backingfile)

test_that("Error: already exists", {
  expect_error(snp_readBed(bedfile, backingfile = path),
               sprintf("File '%s' already exists.", paste0(path, ".bk")),
               fixed = TRUE)
})

################################################################################

test_that("same sign as PLINK (no switch 0 <-> 2)", {
  skip_on_os("solaris"); skip_if_offline("www.cog-genomics.org")
  plink <- download_plink(verbose = FALSE)
  prefix <- sub_bed(bedfile)
  tmp <- tempfile()
  system(glue::glue("{plink} --bfile {prefix} --assoc --allow-no-sex --out {tmp}"),
         ignore.stdout = TRUE, ignore.stderr = TRUE)

  gwas <- big_univLogReg(G, test$fam$affection - 1L)
  sumstats <- bigreadr::fread2(paste0(tmp, ".assoc"))
  expect_gt(cor(gwas$estim, log(sumstats$OR)), 0.99)
})

################################################################################

test_that("snp_readBed2() works", {

  bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")

  ind.row <- sample(nrow(G), nrow(G) / 2, replace = TRUE)
  ind.col <- sample(ncol(G), ncol(G) / 2, replace = TRUE)

  test2 <- snp_attach(snp_readBed2(bedfile, backingfile = tempfile(),
                                   ind.row, ind.col))
  test3 <- snp_attach(subset(test, ind.row, ind.col))
  expect_equal(test2$genotypes[], test3$genotypes[])
  expect_equal(test2[-1], test3[-1])

  test4 <- snp_attach(snp_readBed2(bedfile, backingfile = tempfile()))
  expect_equal(test4$genotypes[], test$genotypes[])
  expect_equal(test4[-1], test[-1])
})

################################################################################
