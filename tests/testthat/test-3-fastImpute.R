################################################################################

context("FAST_IMPUTE")

################################################################################

test_that("fast imputation (xgboost) works", {

  skip_if_not_installed("xgboost")

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

})

################################################################################

test_that("fast imputation (simple) works", {

  G <- snp_attachExtdata("example-missing.bed")$genotypes
  # cbind(colMeans(G[]), colMeans(G[], na.rm = TRUE))
  expect_error(snp_fastImputeSimple(G, "mean"), "should be one of")

  G2 <- snp_fastImputeSimple(G, "mean0")
  expect_equal(G[c(18, 72), 400], rep(NA_real_, 2))
  expect_equal(G2[c(18, 72), 400], rep(1, 2))

  G3 <- snp_fastImputeSimple(G, "mean2")
  expect_equal(G3[c(18, 72), 400], rep(1.01, 2))

  G4 <- snp_fastImputeSimple(G, "mode")
  expect_equal(G4[c(4, 12), 1], rep(0, 2))
  expect_equal(G4[c(18, 72), 400], rep(1, 2))

  imp_val <- replicate(1000, {
    G5 <- snp_fastImputeSimple(G, "random")
    G5[c(18, 72), 400]
  })
  abs_val <- rbinom(2e5, size = 2, prob = mean(G[, 400], na.rm = TRUE) / 2)
  prop <- table(imp_val) / (table(abs_val) / 100)
  expect_equal(as.vector(prop), rep(1, 3), tolerance = 0.1)

})

################################################################################
