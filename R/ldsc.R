################################################################################

# heteroscedasticity and overcounting weights
WEIGHTS <- function(pred, w_ld) {
  1 / (pred^2 * w_ld)
}

crossprod2 <- function(x, y) drop(base::crossprod(x, y))

# equivalent to stats::lm.wfit(cbind(1, x), y, w)
wlm <- function(x, y, w) {
  wx <- w * x
  W   <- sum(w)
  WX  <- sum(wx)
  WY  <- crossprod2(w,  y)
  WXX <- crossprod2(wx, x)
  WXY <- crossprod2(wx, y)
  alpha <- (WXX * WY - WX * WXY) / (W * WXX - WX^2)
  beta  <- (WXY * W  - WX * WY)  / (W * WXX - WX^2)
  list(intercept = alpha, slope = beta, pred = x * beta + alpha)
}

# equivalent to stats::lm.wfit(as.matrix(x), y, w)
wlm_no_int <- function(x, y, w) {
  wx <- w * x
  WXX <- crossprod2(wx, x)
  WXY <- crossprod2(wx, y)
  beta  <- WXY / WXX
  list(slope = beta, pred = x * beta)
}

################################################################################

#' LD score regression
#'
#' @param ld_score Vector of LD scores.
#' @param ld_size Number of variants used to compute `ld_score`.
#' @param chi2 Vector of chi-squared statistics.
#' @param sample_size Sample size of GWAS corresponding to chi-squared statistics.
#'   Possibly a vector, or just a single value.
#' @param chi2_thr1 Threshold on `chi2` in step 1. Default is `30`.
#'   This is equivalent to parameter `--two-step`.
#' @param chi2_thr2 Threshold on `chi2` in step 2. Default is `Inf` (none).
#' @param blocks Either a single number specifying the number of blocks,
#'   or a vector of integers specifying the block number of each `chi2` value.
#'   Default is `200` for `snp_ldsc()`, dividing into 200 blocks of approximately
#'   equal size. `NULL` can also be used to skip estimating standard errors,
#'   which is the default for `snp_ldsc2()`.
#' @param intercept You can constrain the intercept to some value (e.g. 1).
#'   Default is `NULL` in `snp_ldsc()` (the intercept is estimated)
#'   and is `1` in `snp_ldsc2()` (the intercept is fixed to 1).
#'   This is equivalent to parameter `--intercept-h2`.
#' @inheritParams bigsnpr-package
#'
#' @return Vector of 4 values (or only the first 2 if `blocks = NULL`):
#'  - `[["int"]]`: LDSC regression intercept,
#'  - `[["int_se"]]`: SE of this intercept,
#'  - `[["h2"]]`: LDSC regression estimate of (SNP) heritability (also see
#'    [coef_to_liab]),
#'  - `[["h2_se"]]`: SE of this heritability estimate.
#'
#' @importFrom bigassertr assert_one_int
#'
#' @export
#'
snp_ldsc <- function(ld_score, ld_size, chi2, sample_size,
                     blocks = 200,
                     intercept = NULL,
                     chi2_thr1 = 30,
                     chi2_thr2 = Inf,
                     ncores = 1) {

  chi2 <- chi2 + 1e-8
  assert_pos(chi2, strict = TRUE)
  assert_lengths(chi2, ld_score)
  assert_one_int(ld_size)

  M <- length(chi2)
  if (length(sample_size) == 1) {
    sample_size <- rep(sample_size, M)
  } else {
    assert_lengths(sample_size, chi2)
  }

  if (is.null(blocks)) {

    #### step 1 ####

    step1_int <- if (is.null(intercept)) {

      ind_sub1 <- which(chi2 < chi2_thr1)
      w_ld <- pmax(ld_score[ind_sub1], 1)
      x1 <- (ld_score / ld_size * sample_size)[ind_sub1]
      y1 <- chi2[ind_sub1]

      pred0 <- y1
      for (i in 1:100) {
        pred <- wlm(x1, y1, w = WEIGHTS(pred0, w_ld))$pred
        if (max(abs(pred - pred0)) < 1e-6) break
        pred0 <- pred
      }
      wlm(x1, y1, w = WEIGHTS(pred0, w_ld))$intercept

    } else intercept

    #### step 2 ####

    ind_sub2 <- which(chi2 < chi2_thr2)
    w_ld <- pmax(ld_score[ind_sub2], 1)
    x <- (ld_score / ld_size * sample_size)[ind_sub2]
    y <- chi2[ind_sub2]
    yp <- y - step1_int

    pred0 <- y
    for (i in 1:100) {
      pred <- step1_int + wlm_no_int(x, yp, w = WEIGHTS(pred0, w_ld))$pred
      if (max(abs(pred - pred0)) < 1e-6) break
      pred0 <- pred
    }
    step2_h2 <- wlm_no_int(x, yp, w = WEIGHTS(pred0, w_ld))$slope

    c(int = step1_int, h2 = step2_h2)

  } else {

    #### delete-a-group jackknife variance estimator ####

    if (length(blocks) == 1) {
      blocks <- sort(rep_len(seq_len(blocks), M))
    } else {
      assert_lengths(blocks, chi2)
    }
    ind_blocks <- split(seq_along(blocks), blocks)
    h_blocks <- M / lengths(ind_blocks)

    bigparallelr::register_parallel(ncores)

    delete_values <- foreach(ind_rm = c(list(NULL), ind_blocks), .combine = "cbind") %dopar% {
      keep <- which(!seq_along(ld_score) %in% ind_rm)
      snp_ldsc(ld_score[keep], ld_size, chi2[keep], sample_size[keep],
               NULL, intercept, chi2_thr1, chi2_thr2)
    }
    estim <- delete_values[, 1]

    # https://doi.org/10.1023/A:1008800423698
    int_pseudovalues <- h_blocks * estim[1] - (h_blocks - 1) * delete_values[1, -1]
    h2_pseudovalues  <- h_blocks * estim[2] - (h_blocks - 1) * delete_values[2, -1]

    int_J <- sum(int_pseudovalues / h_blocks)
    h2_J  <- sum( h2_pseudovalues / h_blocks)

    c(int    = int_J,
      int_se = sqrt(mean((int_pseudovalues - int_J)^2 / (h_blocks - 1))),
      h2     = h2_J,
      h2_se  = sqrt(mean(( h2_pseudovalues -  h2_J)^2 / (h_blocks - 1))))

  }
}

################################################################################

#' @rdname snp_ldsc
#'
#' @param corr Sparse correlation matrix. Can also be an [SFBM][SFBM-class].
#' @param df_beta A data frame with 3 columns:
#'   - `$beta`: effect size estimates
#'   - `$beta_se`: standard errors of effect size estimates
#'   - `$n_eff`: either GWAS sample size(s) when estimating `beta` for a
#'     continuous trait, or in the case of a binary trait, this is
#'     `4 / (1 / n_control + 1 / n_case)`; in the case of a meta-analysis, you
#'     should sum the effective sample sizes of each study instead of using the
#'     total numbers of cases and controls, see \doi{10.1016/j.biopsych.2022.05.029};
#'     when using a mixed model, the effective sample size needs to be adjusted
#'     as well, see \doi{10.1016/j.xhgg.2022.100136}.
#' @param ind.beta Indices in `corr` corresponding to `df_beta`. Default is all.
#' @inheritDotParams snp_ldsc chi2_thr1 chi2_thr2
#'
#' @export
#'
#' @examples
#' bigsnp <- snp_attachExtdata()
#' G <- bigsnp$genotypes
#' y <- bigsnp$fam$affection - 1
#' corr <- snp_cor(G, ind.col = 1:1000)
#'
#' gwas <- big_univLogReg(G, y, ind.col = 1:1000)
#' df_beta <- data.frame(beta = gwas$estim, beta_se = gwas$std.err,
#'                       n_eff = 4 / (1 / sum(y == 0) + 1 / sum(y == 1)))
#'
#' snp_ldsc2(corr, df_beta)
#' snp_ldsc2(corr, df_beta, blocks = 20, intercept = NULL)
#'
snp_ldsc2 <- function(corr, df_beta,
                      blocks = NULL,
                      intercept = 1,
                      ncores = 1,
                      ind.beta = cols_along(corr),
                      ...) {

  assert_df_with_names(df_beta, c("beta", "beta_se", "n_eff"))
  assert_lengths(ind.beta, rows_along(df_beta))
  stopifnot(all(ind.beta %in% cols_along(corr)))
  assert_pos(df_beta$beta_se, strict = TRUE)

  full_ld <- if (inherits(corr, "dsCMatrix")) {
    sp_colSumsSq_sym(corr@p, corr@i, corr@x)
  } else if (inherits(corr, "SFBM")) {
    ld_scores_sfbm(corr, ind_sub = cols_along(corr) - 1L, ncores = ncores)
  } else {
    Matrix::colSums(corr^2)
  }

  snp_ldsc(
    ld_score    = full_ld[ind.beta],
    ld_size     = ncol(corr),
    chi2        = (df_beta$beta / df_beta$beta_se)^2,
    sample_size = df_beta$n_eff,
    blocks      = blocks,
    intercept   = intercept,
    ncores      = ncores,
    ...
  )
}

################################################################################

#' Liability scale
#'
#' Coefficient to convert to the liability scale. E.g. h2_liab = coef * h2_obs.
#'
#' @param K_pop Prevalence in the population.
#' @param K_gwas Prevalence in the GWAS. You should provide this if you used
#'   (`n_case + n_control`) as sample size. If using the effective sample size
#'   `4 / (1 / n_case + 1 / n_control)` instead, you should keep the default
#'   value of `K_gwas = 0.5` as the GWAS case-control ascertainment is already
#'   accounted for in the effective sample size.
#'
#' @return Scaling coefficient to convert e.g. heritability to the liability scale.
#' @export
#'
#' @examples
#' h2 <- 0.2
#' h2 * coef_to_liab(0.02)
coef_to_liab <- function(K_pop, K_gwas = 0.5) {

  z <- stats::dnorm(stats::qnorm(min(K_pop, 1 - K_pop)))

  (K_pop * (1 - K_pop) / z)^2 / (K_gwas * (1 - K_gwas))
}

################################################################################
