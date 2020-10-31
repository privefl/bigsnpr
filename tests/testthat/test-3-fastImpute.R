################################################################################

context("FAST_IMPUTE")

options(bigstatsr.check.parallel.blas = FALSE)

################################################################################

test_that("fast imputation (xgboost) works", {

  skip_if_not_installed("xgboost")
  skip_if(is_cran)

  suppressMessages({
    library(Matrix)
    library(foreach)
    library(magrittr)
    library(dplyr)
  })

  bigsnp <- snp_attachExtdata()
  G <- bigsnp$genotypes
  expect_equal(G$code256, CODE_012)
  corr <- snp_cor(G)
  ind <- which(apply(corr, 2, function(x) max(x[x < 1])) > 0.6)

  indNA <- foreach(j = ind, .combine = "rbind") %do% {
    cbind(sample(nrow(G), 100), j)
  }

  # Create a copy of our bigSNP example
  tmpfile <- tempfile()
  bigsnp.copy <- snp_writeBed(bigsnp, bedfile = paste0(tmpfile, ".bed")) %>%
    snp_readBed(backingfile = tmpfile) %>%
    snp_attach()
  # Fill some missing values
  GNA <- bigsnp.copy$genotypes
  GNA[indNA] <- as.raw(3)
  counts <- big_counts(GNA)
  expect_equal(sum(counts[4, ]), nrow(indNA))

  ################################################################################

  # Fast imputation
  time <- system.time(
    infos <- snp_fastImpute(GNA, infos.chr = bigsnp$map$chromosome, p.train = 0.6)
  )
  # expect_lt(time[3], 10)

  # Still NAs
  counts <- big_counts(GNA)
  expect_equal(sum(counts[4, ]), nrow(indNA))

  GNA$code256 <- CODE_IMPUTE_PRED
  infosNA <- tibble(
    col = indNA[, 2],
    error = (GNA[indNA] != G[indNA])
  ) %>%
    group_by(col) %>%
    summarise(nb_err = sum(error)) %>%
    mutate(nb_err_est = (infos[1, ] * infos[2, ] * nrow(G))[col])

  pval <- anova(lm(nb_err ~ nb_err_est - 1, data = infosNA))$`Pr(>F)`[[1]]
  expect_lt(pval, 2e-16)

  ################################################################################

  fake <- snp_fake(100, 200)
  G <- fake$genotypes
  G[] <- sample(as.raw(0:3), size = length(G), replace = TRUE)
  CHR <- fake$map$chromosome
  # 'infos.chr' needs to have the same length as the number of SNPs
  expect_error(snp_fastImpute(G, 1), "Incompatibility between dimensions.",
               fixed = TRUE)
  # you can't impute randomness, right?
  expect_gt(mean(snp_fastImpute(G, CHR)[2, ]), 0.5)

  ################################################################################

  fake <- snp_fake(100, 200)
  G <- fake$genotypes
  G[] <- sample(as.raw(0:3), size = length(G), replace = TRUE)
  CHR <- sort(sample(1:3, ncol(G), TRUE))
  G2 <- bigsnpr:::FBM_infos(G)
  ind <- 1:11 + 150
  G2[1, ind] <- seq(0, 1, 0.1)
  G2 <- snp_fastImpute(G, CHR)
  G3 <- G$copy(CODE_IMPUTE_PRED)
  nbNA <- colSums(is.na(G3[]))
  expect_true(all(nbNA[ind] > 0))
  expect_true(all(nbNA[-ind] == 0))

  ################################################################################

  bigsnp <- snp_attachExtdata()
  G <- bigsnp$genotypes
  CHR <- bigsnp$map$chromosome

  elemNA <- sample(length(G), size = 20)
  G2 <- big_copy(G); G2[elemNA] <- as.raw(3)  # NA
  G3 <- big_copy(G); G3[elemNA] <- as.raw(3)  # NA

  expect_equal(snp_fastImpute(G2, CHR, seed = 2)[],
               snp_fastImpute(G3, CHR, seed = 2, ncores = 2)[])
  expect_equal(G2[], G3[])

})

################################################################################

test_that("fast imputation (simple) works", {

  G <- snp_attachExtdata("example-missing.bed")$genotypes
  # cbind(colMeans(G[]), colMeans(G[], na.rm = TRUE))
  expect_error(snp_fastImputeSimple(G, "mean"), "should be one of")

  expect_warning(G1 <- snp_fastImputeSimple(G, "zero"), "deprecated")
  expect_equal(G[c(18, 72), 400], rep(NA_real_, 2))
  expect_equal(G1[c(18, 72), 400], rep(0, 2))

  G2 <- snp_fastImputeSimple(G, "mean0")
  expect_equal(G[c(18, 72), 400], rep(NA_real_, 2))
  expect_equal(G2[c(18, 72), 400], rep(1, 2))

  G3 <- snp_fastImputeSimple(G, "mean2")
  expect_equal(G3[c(18, 72), 400], rep(1.01, 2))

  G4 <- snp_fastImputeSimple(G, "mode")
  expect_equal(G4[c(4, 12), 1], rep(0, 2))
  expect_equal(G4[c(18, 72), 400], rep(1, 2))

  imp_val <- replicate(500, {
    G5 <- snp_fastImputeSimple(G, "random")
    G5[c(18, 72), 400]
  })
  p <- mean(G[, 400], na.rm = TRUE) / 2
  prob <- c((1 - p)^2, 2 * p * (1 - p), p^2)
  expect_gt(stats::chisq.test(table(imp_val), p = prob)$p.value, 1e-4)

  expect_equal(snp_fastImputeSimple(G, method = "mean2")[],
               snp_fastImputeSimple(G, method = "mean2", ncores = 2)[])
})

################################################################################
