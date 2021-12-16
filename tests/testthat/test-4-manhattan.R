context("MANHATTAN")

test_that("snp_manhattan() works with unordered data", {

  test <- snp_attachExtdata()
  G <- test$genotypes

  gwas <- big_univLogReg(G, test$fam$affection - 1L)
  N <- ncol(G)
  CHR <- sort(rep_len(1:2, N))
  POS <- 1:N * 1000
  rand <- sample(N)
  # plot_grid(
  #   snp_manhattan(gwas, CHR, POS, dist.sep.chrs = 0),
  #   snp_manhattan(gwas[rand, ], CHR[rand], POS[rand], dist.sep.chrs = 0)
  # )

  expect_equal(snp_manhattan(gwas, CHR, POS),
               snp_manhattan(gwas[rand, ], CHR[rand], POS[rand]))
  expect_equal(snp_manhattan(gwas, CHR, POS, npoints = 500),
               snp_manhattan(gwas[rand, ], CHR[rand], POS[rand], npoints = 500))
  expect_equal(snp_manhattan(gwas, as.character(CHR), POS, npoints = 500),
               snp_manhattan(gwas[rand, ], as.character(CHR[rand]), POS[rand], npoints = 500))
})
