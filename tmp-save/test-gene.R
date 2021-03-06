################################################################################

context("GENE")

################################################################################

rsid <- c("rs3934834", "rs3737728", "rs6687776", "rs9651273", "rs4970405",
          "rs12726255", "rs2298217", "rs4970362", "rs9660710", "rs4970420")
genes <- c("LOC105378948:105378948,RNF223:401934", rep("C1orf159:54991", 5),
           rep(NA, 4))

test_that("Getting genes works", {
  expect_equal(snp_gene(rsid), genes)
  expect_equal(snp_gene(rsid, ncores = 2), genes)
})

################################################################################
