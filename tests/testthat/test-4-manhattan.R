context("MANHATTAN")

test_that("snp_manhattan() works with unordered data", {

  test <- snp_attachExtdata()
  G <- test$genotypes

  gwas <- big_univLogReg(G, test$fam$affection - 1L)
  N <- ncol(G)
  CHR <- sort(rep_len(1:2, N))
  POS <- 1:N * 1000
  expect_error(snp_manhattan(gwas, CHR[1:5], POS[1:5]),
               "Incompatibility between dimensions")

  rand <- sample(N)
  # plot_grid(
  #   snp_manhattan(gwas, CHR, POS, dist.sep.chrs = 0),
  #   snp_manhattan(gwas[rand, ], CHR[rand], POS[rand], dist.sep.chrs = 0)
  # )

  expect_same_plot <- function(p1, p2) {
    png1 <- ggplot2::ggsave(tempfile(fileext = ".png"), p1, width = 8, height = 6)
    png2 <- ggplot2::ggsave(tempfile(fileext = ".png"), p2, width = 8, height = 6)
    expect_identical(readBin(png1, what = raw(), n = 1e6),
                     readBin(png2, what = raw(), n = 1e6))
  }

  expect_failure(expect_same_plot(snp_manhattan(gwas, CHR, POS),
                                  snp_manhattan(gwas, CHR, POS, npoints = 500)))

  expect_same_plot(snp_manhattan(gwas, CHR, POS),
                   snp_manhattan(gwas[rand, ], CHR[rand], POS[rand]))
  expect_same_plot(snp_manhattan(gwas, CHR, POS, npoints = 500),
                   snp_manhattan(gwas[rand, ], CHR[rand], POS[rand], npoints = 500))
  expect_same_plot(snp_manhattan(gwas, as.character(CHR), POS, npoints = 500),
                   snp_manhattan(gwas[rand, ], as.character(CHR[rand]), POS[rand], npoints = 500))
})
