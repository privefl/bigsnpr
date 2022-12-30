################################################################################

context("LASSOSUM2")

################################################################################

test_that("lassosum2 works", {

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
  df_beta <- snp_match(sumstats, map, join_by_pos = FALSE)[c("beta", "beta_se",
                                                             "n_eff", "_NUM_ID_")]

  ind_var <- df_beta$`_NUM_ID_`
  corr0 <- snp_cor(G, ind.col = ind_var, size = 3 / 1000,
                   infos.pos = POS2[ind_var], ncores = 2)
  corr <- as_SFBM(corr0, compact = sample(c(TRUE, FALSE), 1))
  rm(corr0)

  # lassosum2
  nlam <- sample(8:15, 1)
  beta_grid <- snp_lassosum2(corr, df_beta, nlambda = nlam, maxiter = 100, ncores = 2)
  pred_grid <- big_prodMat(G, beta_grid, ind.col = ind_var)
  expect_gt(max(cor(pred_grid, y), na.rm = TRUE), 0.4)
  params <- attr(beta_grid, "grid_param")
  expect_equal(nrow(params), 4 * nlam)
  expect_true(all(params$num_iter <= 101))

  beta_grid2 <- snp_lassosum2(corr, df_beta, nlambda = nlam, maxiter = 100, ncores = 2)
  attr(beta_grid2, "grid_param")$time <- attr(beta_grid, "grid_param")$time <- NULL
  expect_identical(beta_grid2, beta_grid)  # no sampling (deterministic), so reproducible

  # ind.corr
  ind.sub <- sample(ncol(corr), 30e3)
  beta_grid2 <- snp_lassosum2(as_SFBM(corr[ind.sub, ind.sub]), df_beta[ind.sub, ],
                              nlambda = nlam, maxiter = 50, ncores = 2)
  beta_grid3 <- snp_lassosum2(corr, df_beta[ind.sub, ], ind.corr = ind.sub,
                              nlambda = nlam, maxiter = 50, ncores = 2)
  expect_equal(beta_grid2, beta_grid3, check.attributes = FALSE)
  # plot(attr(beta_grid2, "grid_param")$time, attr(beta_grid3, "grid_param")$time)
  # lm(attr(beta_grid3, "grid_param")$time ~ attr(beta_grid2, "grid_param")$time + 0)$coef[[1]]

  # Errors
  expect_error(snp_lassosum2(corr, df_beta[-1]),
               "'df_beta' should have element 'beta'.")
  expect_error(snp_lassosum2(corr, df_beta[-2]),
               "'df_beta' should have element 'beta_se'.")
})

################################################################################
