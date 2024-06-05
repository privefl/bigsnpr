################################################################################

#' Quality control and imputation of summary statistics
#'
#' Quality control of summary statistics, by comparing z-scores with z-scores
#' imputed using nearby correlated variants. Also performs imputation.
#'
#' @inheritParams snp_ldpred2_auto
#' @param sd0 Vector of standard deviations (of genotypes). You can estimate
#'   those using `sqrt(2 * af * (1 - af) * INFO)`, where `af` are the allele
#'   frequencies and `INFO` the imputation info scores.
#' @param thr_highld Threshold on squared correlations to decide which variants
#'   to use for imputation. Default is `0.2`. This MUST be equal (or larger)
#'   than `thr_r2` used in `snp_cor_extendedThr()` to compute `corr`.
#' @param thr_ld Threshold on LD scores for restricting variants in LD to use
#'   for imputation. This is used to avoid redundancy and increase speed.
#'   Default is `30`.
#' @param removed Logical vector to discard some variants from being used in
#'   the imputation of other variants. However, these variants will be imputed.
#' @param max_run Maximum number of iterations (and therefore the maximum number
#'   of variants that can be removed). Default is the number of input variants.
#' @param pval_thr P-value threshold for QCing. Default is `5e-8`.
#' @param print_iter Whether to print the iteration numbers. Default is `FALSE`.
#'
#' @import foreach
#'
#' @return A tibble (data frame) with 8 variables for each variants:
#'   - `$beta_imp`: imputed effect sizes,
#'   - `$beta_se_imp`: imputed standard errors (corresponding to `$n_eff_imp`),
#'   - `$n_eff_imp`: imputed sample sizes,
#'   - `$r2_imp`: multiple $R^2$ from imputation,
#'   - `$r2_max`: maximum $r^2$ with variants used for imputing,
#'   - `$chi2_qc`: the chi-squared statistics,
#'   - `$pval_qc`: the corresponding p-values (df = 1),
#'   - `$rm_qc`: whether the p-values are significant (or `removed`).
#' @export
#'
snp_qcimp_sumstats <- function(corr, df_beta, sd0,
                               thr_highld = 0.2,
                               thr_ld = 30,
                               removed = rep(FALSE, nrow(df_beta)),
                               pval_thr = 5e-8,
                               max_run = nrow(df_beta),
                               print_iter = FALSE,
                               ncores = 1) {

  assert_df_with_names(df_beta, c("beta", "beta_se", "n_eff"))
  df_beta <- df_beta[, c("beta", "beta_se", "n_eff")]

  assert_lengths(rows_along(corr), cols_along(corr),
                 rows_along(df_beta), sd0, removed)
  assert_nona(df_beta)
  assert_nona(sd0)

  keep <- !removed | is.na(removed)

  sd <- with(df_beta, 1 / sqrt(n_eff * beta_se^2 + beta^2))
  opt_scaling <- stats::optimize(
    function(scaling) sum(abs(sd0 - sd * scaling)[keep]),
    interval = quantile((sd0 / sd)[keep], c(0.001, 0.999))
  )$min
  sd <- sd * opt_scaling

  diag_scaler <- with(df_beta, max(n_eff) / n_eff)

  chi2_thr <- stats::qchisq(pval_thr, df = 1, lower.tail = FALSE)
  thr <- rep(thr_highld, nrow(df_beta))

  # preallocate results
  all_beta     <- rep(NA_real_,   nrow(df_beta))
  all_neff     <- rep(NA_real_,   nrow(df_beta))
  all_max_r2   <- rep(NA_real_,   nrow(df_beta))
  all_multi_r2 <- rep(NA_real_,   nrow(df_beta))
  all_chi2     <- rep(-1,         nrow(df_beta))
  all_id_high  <- rep(list(NULL), nrow(df_beta))

  keep_high <- rep(FALSE, nrow(df_beta))  # placeholder

  # for the first iteration, impute all
  id_this <- rows_along(df_beta)

  bigparallelr::register_parallel(ncores)

  for (i_run in seq_len(max_run)) {

    if (print_iter) cat(i_run, "")

    if (length(id_this) > 0) {

      thr[] <- thr_highld
      thr[which(all_chi2 > (2 * chi2_thr))] <- Inf

      FUNs <- c("find_ld_friends", "test_ld_scores", "cor2cov_inplace")
      all_res_this <- foreach(id_current = id_this, .export = FUNs) %dopar% {

        # make sure use local copies before C++
        keep[1] <- keep[1]
        keep_high[1] <- FALSE

        # all LD friends of current
        high_ld <- find_ld_friends(corr, id_current - 1L, keep, thr)
        # variant selection, to avoid redundancy and to fasten
        test_ld_scores(corr, ord = order(high_ld[[2]]^2, decreasing = TRUE) - 1L,
                       ind = high_ld[[1]], keep = keep_high, thr = thr_ld)
        id_high <- which(keep_high)
        keep_high[id_high] <- FALSE  # reset to all FALSE
        if (length(id_high) == 0)
          return(list(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NULL))

        sd_high    <- sd[id_high]
        sd_current <- sd[id_current]
        cor_high_current <- corr$dense_acc(id_high, id_current)
        cov_high_current <- cor_high_current * sd_high * sd_current
        cov_high_high <- cor2cov_inplace(corr$dense_acc(id_high, id_high),
                                         sd = sd_high)

        reg <- if (i_run == 1) 0.5 else 0.001
        diag(cov_high_high) <- diag(cov_high_high) * diag_scaler[id_high] + reg
        w <- solve(cov_high_high, cov_high_current)

        multi_R2 <- crossprod(cov_high_current, w) / sd_current^2
        max_R2 <- max(cor_high_current^2)
        R2 <- min(multi_R2, max_R2)

        n_eff_r2 <- df_beta$n_eff[id_high] * cor_high_current**2
        n_eff_imp <- min(stats::weighted.mean(n_eff_r2, abs(w)), max(n_eff_r2))

        beta_imp <- crossprod(df_beta$beta[id_high] * sd_high^2, w) / sd_current^2

        diff_z <- (beta_imp - df_beta$beta[id_current]) / df_beta$beta_se[id_current]
        z_scaler <- max(with(df_beta[id_high, ], abs(beta / beta_se))) / 4
        chi2 <- (diff_z / max(z_scaler, 1))^2 / max(1 - R2, 0.01) * R2

        list(beta_imp, n_eff_imp, max_R2, multi_R2, chi2, id_high)
      }

      # save current results
      all_beta[id_this]     <- sapply(all_res_this, function(x) x[[1]])
      all_neff[id_this]     <- sapply(all_res_this, function(x) x[[2]])
      all_max_r2[id_this]   <- sapply(all_res_this, function(x) x[[3]])
      all_multi_r2[id_this] <- sapply(all_res_this, function(x) x[[4]])
      all_chi2[id_this]     <- sapply(all_res_this, function(x) x[[5]])
      all_id_high[id_this]  <- lapply(all_res_this, function(x) x[[6]])
    }

    if (i_run > 1) {  # redo them all in the 2nd iter, for updating `thr`
      # get the worst variant and remove it
      id_worst <- which.max(all_chi2 * keep)
      if (all_chi2[id_worst] > chi2_thr) {
        keep[id_worst] <- FALSE
        # find which variants used id_worst for imputation, and update them
        run_next <- sapply(all_id_high, function(id_high) id_worst %in% id_high)
        id_this <- which(run_next)
      } else break
    }

  }

  pval <- stats::pchisq(all_chi2, df = 1, lower.tail = FALSE)
  # note that a variant that has been removed in one iteration can still have
  # a non-significant p-value in the end; this can prevent some false positives

  tibble::tibble(
    beta_imp    = all_beta,
    beta_se_imp = sqrt((opt_scaling^2 / sd0^2 - all_beta^2) / all_neff),
    n_eff_imp   = all_neff,
    r2_imp      = all_multi_r2,
    r2_max      = all_max_r2,
    chi2_qc     = all_chi2,
    pval_qc     = pval,
    rm_qc       = removed | (pval < pval_thr)
  )
}

################################################################################
