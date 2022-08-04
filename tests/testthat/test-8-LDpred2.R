################################################################################

context("LDPRED2")

################################################################################

test_that("LDpred2 works", {

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
  corr0 <- snp_cor(G, ind.col = ind_var, size = 3 / 1000,
                   infos.pos = POS2[ind_var], ncores = 2)
  ld <- bigsnpr:::sp_colSumsSq_sym(p = corr0@p, i = corr0@i, x = corr0@x)
  corr <- as_SFBM(corr0, compact = sample(c(TRUE, FALSE), 1))
  rm(corr0)

  # LD score regression
  (ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                  sample_size = n_eff, blocks = NULL)))

  # LDpred2-inf
  beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = ldsc[["h2"]])
  pred_inf <- big_prodVec(G, beta_inf, ind.col = df_beta[["_NUM_ID_"]])
  expect_gt(cor(pred_inf, y), 0.2)

  # LDpred2-gibbs
  p_seq <- signif(seq_log(1e-3, 1, length.out = 7), 1)
  params <- expand.grid(p = p_seq, h2 = ldsc[["h2"]], sparse = c(FALSE, TRUE))
  expect_equal(dim(params), c(14, 3))
  expect_equal(dim(snp_ldpred2_grid(corr, df_beta, params[1, ])), c(nrow(df_beta), 1))
  beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = 2)
  pred_grid <- big_prodMat(G, beta_grid, ind.col = df_beta[["_NUM_ID_"]])
  expect_gt(max(cor(pred_grid, y)), 0.4)

  # LDpred2-uncertainty
  expect_error(snp_ldpred2_grid(corr, df_beta, params, return_sampling_betas = TRUE),
               "Only one set of parameters is allowed")
  beta_sample <- snp_ldpred2_grid(corr, df_beta, params[3, ], num_iter = 200,
                                  return_sampling_betas = TRUE)
  expect_equal(dim(beta_sample), c(nrow(df_beta), 200))
  expect_gt(cor(rowMeans(beta_sample), beta_grid[, 3]), 0.9)

  # LDpred2-auto
  beta_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = ldsc[["h2"]],
                                burn_in = 200, num_iter = 200, sparse = TRUE,
                                shrink_corr = runif(1, 0.9, 1))
  expect_length(beta_auto, 1)
  mod <- beta_auto[[1]]
  expect_null(dim(mod$beta_est))
  pred_auto <- big_prodVec(G, mod$beta_est, ind.col = df_beta[["_NUM_ID_"]])

  expect_gt(cor(pred_auto, y), 0.4)
  expect_gt(cor(mod$beta_est_sparse, mod$beta_est), 0.9)
  expect_lt(mean(mod$beta_est == 0), 0.001)
  expect_gt(mean(mod$beta_est_sparse == 0), 0.5)
  expect_equal(mod$h2_init, ldsc[["h2"]])
  expect_equal(mod$p_init, 0.1)
  expect_equal(mod$h2_est, mean(tail(mod$path_h2_est, 200)))
  expect_equal(mod$p_est,  mean(tail(mod$path_p_est,  200)))
  expect_length(mod$postp_est, length(mod$beta_est))
  expect_null(dim(mod$postp_est))
  expect_equal(mean(mod$postp_est), mod$p_est, tolerance = 0.01)
  expect_equal(dim(mod$sample_beta), c(ncol(corr), 0))
  beta_hat <- with(df_beta, beta / sqrt(n_eff * beta_se^2 + beta^2))
  expect_gt(cor(mod$corr_est, beta_hat), 0.3)

  beta_auto2 <- snp_ldpred2_auto(corr, df_beta, h2_init = ldsc[["h2"]],
                                 burn_in = 200, num_iter = 200,
                                 vec_p_init = 1, allow_jump_sign = FALSE)
  expect_gt(cor(beta_auto2[[1]]$beta_est, mod$beta_est), 0.9)
  expect_lt(beta_auto2[[1]]$p_est, 0.1)

  # Sampling betas
  beta_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = ldsc[["h2"]],
                                burn_in = 200, num_iter = 200, report_step = 10)
  bsamp <- beta_auto[[1]]$sample_beta
  expect_equal(dim(bsamp), c(ncol(corr), 20))
  expect_s4_class(bsamp, "dgCMatrix")
  h2_est <- apply(bsamp, 2, function(x) crossprod(x, bigsparser::sp_prodVec(corr, x)))
  expect_equal(h2_est, beta_auto[[1]]$path_h2_est[seq(210, 400, by = 10)])

  # Errors
  expect_error(snp_ldpred2_inf(corr, df_beta[-6]),
               "'df_beta' should have element 'beta'.")
  expect_error(snp_ldpred2_grid(corr, df_beta[-6], params),
               "'df_beta' should have element 'beta'.")
  expect_error(snp_ldpred2_auto(corr, df_beta[-6]),
               "'df_beta' should have element 'beta'.")

  # Parallelism
  all_beta <- replicate(
    n = 10, simplify = FALSE,
    snp_ldpred2_grid(corr, df_beta, params[c(1, 8), ], ncores = 2))
  all_beta <- replicate(
    n = 5, simplify = FALSE,
    snp_ldpred2_auto(corr, df_beta, h2_init = ldsc[["h2"]],
                     vec_p_init = seq_log(1e-4, 0.9, 6),
                     burn_in = 100, num_iter = 100,
                     sparse = TRUE, ncores = 2))

  # Reproducibility

  set.seed(1)
  grid1 <- snp_ldpred2_grid(corr, df_beta, params, ncores = 2)
  grid2 <- snp_ldpred2_grid(corr, df_beta, params, ncores = 2)
  expect_false(identical(grid2, grid1))
  set.seed(1)
  grid3 <- snp_ldpred2_grid(corr, df_beta, params, ncores = 2)
  expect_identical(grid3, grid1)

  set.seed(1)
  beta1 <- snp_ldpred2_grid(corr, df_beta, params[3, ], return_sampling_betas = TRUE)
  set.seed(1)
  beta2 <- snp_ldpred2_grid(corr, df_beta, params[3, ], return_sampling_betas = TRUE)
  expect_identical(beta2, beta1)

  set.seed(1)
  auto1 <- snp_ldpred2_auto(corr, df_beta,
                            burn_in = 50, num_iter = 100, sparse = TRUE,
                            h2_init = 0.3, vec_p_init = seq_log(1e-5, 1, 8),
                            report_step = 5, ncores = 2)
  set.seed(1)
  auto2 <- snp_ldpred2_auto(corr, df_beta,
                            burn_in = 50, num_iter = 100, sparse = TRUE,
                            h2_init = 0.3, vec_p_init = seq_log(1e-5, 1, 8),
                            report_step = 5, ncores = 2)
  expect_identical(auto2, auto1)

  inf1 <- snp_ldpred2_inf(corr, df_beta, h2 = 0.3)
  inf2 <- snp_ldpred2_inf(corr, df_beta, h2 = 0.3)
  expect_identical(inf2, inf1)  # no sampling, so reproducible
})

################################################################################
