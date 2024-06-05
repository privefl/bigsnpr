################################################################################

context("QCIMP_SUMSTATS")

################################################################################

test_that("snp_qcimp_sumstats() works", {

  skip_if(is_cran)
  skip_if_offline("raw.githubusercontent.com")

  bedfile <- file.path(tempdir(), "tmp-data/public-data3.bed")
  if (!file.exists(rdsfile <- sub_bed(bedfile, ".rds"))) {
    zip <- tempfile(fileext = ".zip")
    download.file(
      "https://github.com/privefl/bigsnpr/blob/master/data-raw/public-data3.zip?raw=true",
      destfile = zip, mode = "wb")
    unzip(zip, exdir = tempdir())
    rds <- snp_readBed(bedfile)
    expect_identical(normalizePath(rds), normalizePath(rdsfile))
  }

  obj.bigSNP <- snp_attach(rdsfile)
  G <- obj.bigSNP$genotypes
  y <- obj.bigSNP$fam$affection
  POS2 <- obj.bigSNP$map$genetic.dist + 1000 * obj.bigSNP$map$chromosome

  sumstats <- bigreadr::fread2(file.path(tempdir(), "tmp-data/public-data3-sumstats.txt"))
  sumstats$n_eff <- sumstats$N
  map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
  df_beta <- snp_match(sumstats, map, join_by_pos = FALSE)

  ind_var <- df_beta$`_NUM_ID_`
  corr0 <- snp_cor_extendedThr(
    G, ind.col = ind_var, thr_r2 = 0.2,
    infos.pos = POS2[ind_var], size = 1 / 1000, ncores = 2)
  corr <- as_SFBM(corr0)

  maf <- snp_MAF(G, ind.col = ind_var)
  sd0 <- sqrt(2 * maf * (1 - maf))
  res <- snp_qcimp_sumstats(corr, df_beta, sd0 = sd0, ncores = 2)
  expect_false(any(res$rm_qc, na.rm = TRUE))

  df_beta2 <- df_beta
  ind <- sample(which(!is.na(res$rm_qc) & df_beta$p < 0.001), 20)
  df_beta2$beta[ind] <- -df_beta$beta[ind]
  res2 <- snp_qcimp_sumstats(corr, df_beta2, sd0 = sd0, ncores = 2)

  cor_nona <- function(x, y) cor(x, y, use = "pairwise.complete.obs")

  expect_true(any(res2$rm_qc, na.rm = TRUE))
  expect_lt(median(res2$pval_qc[ind], na.rm = TRUE),
            median(res2$pval_qc, na.rm = TRUE))
  expect_lt(cor_nona(res2$beta_imp[ind], df_beta2$beta[ind]), 0)
  expect_gt(cor_nona(res2$beta_imp, df_beta$beta), 0.5)
  expect_gt(cor_nona(res2$n_eff_imp, res2$r2_imp), 0.9)
  expect_gt(cor_nona(res2$r2_imp, res2$r2_max), 0.9)
  expect_gt(min(res2$beta_se_imp, na.rm = TRUE), 0)
  expect_gt(min(res2$r2_imp, na.rm = TRUE), 0.19)
  expect_gte(min(res2$r2_max, na.rm = TRUE), 0.2)
  expect_lte(max(res2$r2_imp, na.rm = TRUE), 1)
  expect_lte(max(res2$r2_max, na.rm = TRUE), 1)

  # The algorithm should be deterministic
  expect_equal(snp_qcimp_sumstats(corr, df_beta,  sd0 = sd0), res)
  expect_equal(snp_qcimp_sumstats(corr, df_beta2, sd0 = sd0), res2)
})

################################################################################

test_that("cor2cov_inplace() works", {
  x <- matrix(rnorm(100), 10)
  r <- cor(x)
  s <- apply(x, 2, sd)
  r2 <- bigsnpr:::cor2cov_inplace(r, s)
  expect_equal(r2, r)
  expect_equal(diag(r2), s^2)
  expect_equal(r2, cov(x))
})

################################################################################
