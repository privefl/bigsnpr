################################################################################

context("LASSOSUM2")

################################################################################

test_that("lassosum2 works", {

  skip_if(is_cran)

  unzip(test_path("testdata/public-data3.zip"), exdir = tempdir())
  rds <- snp_readBed(file.path(tempdir(), "tmp-data/public-data3.bed"))

  obj.bigSNP <- snp_attach(rds)
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
  ld <- bigsnpr:::sp_colSumsSq_sym(p = corr0@p, i = corr0@i, x = corr0@x)
  corr <- as_SFBM(corr0, compact = sample(c(TRUE, FALSE), 1))
  rm(corr0)

  # lassosum2
  nlam <- sample(8:15, 1)
  beta_grid <- snp_lassosum2(corr, df_beta, nlambda = nlam, ncores = 2)
  pred_grid <- big_prodMat(G, beta_grid, ind.col = ind_var)
  expect_gt(max(cor(pred_grid, y), na.rm = TRUE), 0.4)
  params <- attr(beta_grid, "grid_param")
  expect_equal(nrow(params), 6 * nlam)
  expect_true(all(params$num_iter <= 501))

  beta_grid2 <- snp_lassosum2(corr, df_beta, nlambda = nlam, ncores = 2)
  attr(beta_grid2, "grid_param")$time <- attr(beta_grid, "grid_param")$time <- NULL
  expect_identical(beta_grid2, beta_grid)

  # Errors
  expect_error(snp_lassosum2(corr, df_beta[-1]),
               "'df_beta' should have element 'beta'.")
  expect_error(snp_lassosum2(corr, df_beta[-2]),
               "'df_beta' should have element 'beta_se'.")
})

################################################################################
