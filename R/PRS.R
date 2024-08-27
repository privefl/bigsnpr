################################################################################

prodVecRev <- function(G, betas.col, same.col, ind.row, ind.col) {
  betas.col.mod <- (2 * same.col - 1) * betas.col
  big_prodVec(G, betas.col.mod, ind.row, ind.col) +
    2 * sum(betas.col[!same.col])
}

################################################################################

#' PRS
#'
#' Polygenic Risk Scores with possible clumping and thresholding.
#'
#' @inheritParams bigsnpr-package
#' @param betas.keep Numeric vector of weights associated with each SNP
#'   corresponding to `ind.keep`. You may want to see [bigstatsr::big_univLinReg]
#'   or [bigstatsr::big_univLogReg].
#' @param ind.test The individuals on whom to project the scores. Default uses all.
#' @param ind.keep Column (SNP) indices to use (if using clumping, the
#'   output of [snp_clumping]). Default doesn't clump.
#' @param same.keep A logical vector associated with `betas.keep` whether the
#'   reference allele is the same for G. Default is all `TRUE` (for example when
#'   you train the betas on the same dataset). Otherwise, use [same_ref].
#' @param lpS.keep Numeric vector of `-log10(p-value)` associated with
#'   `betas.keep`. Default doesn't use thresholding.
#' @param thr.list Threshold vector on `lpS.keep` at which SNPs are excluded if
#'   they are not significant enough. Default doesn't use thresholding.
#'
#' @return A matrix of scores, where rows correspond to `ind.test` and
#'   columns correspond to `thr.list`.
#' @export
#'
#' @example examples/example-PRS.R
#'
snp_PRS <- function(G, betas.keep,
                    ind.test = rows_along(G),
                    ind.keep = cols_along(G),
                    same.keep = rep(TRUE, length(ind.keep)),
                    lpS.keep = NULL,
                    thr.list = 0) {

  check_args()
  assert_nona(same.keep)
  assert_type(same.keep, "logical")
  assert_lengths(same.keep, ind.keep)
  assert_lengths(betas.keep, ind.keep)

  # thresholding and projecting
  if (is.null(lpS.keep) || identical(thr.list, 0)) {
    message("'lpS.keep' or 'thr.list' was not specified. Thresholding disabled.")
    scores.all <- as.matrix(prodVecRev(G, betas.keep, same.keep,
                                       ind.test, ind.keep))
  } else {
    assert_lengths(lpS.keep, ind.keep)
    assert_pos(lpS.keep, strict = FALSE)

    scores.all <- matrix(NA_real_, length(ind.test), length(thr.list))
    ind.rem <- seq_along(ind.keep)
    last <- 0
    for (i in order(thr.list, decreasing = TRUE)) {
      pass.thr <- (lpS.keep[ind.rem] > thr.list[i])
      ind <- ind.rem[pass.thr]
      # build score on top of previous ones
      scores.all[, i] <- last <- last +
        prodVecRev(G, betas.keep[ind], same.keep[ind], ind.test, ind.keep[ind])
      # remaining indices amongst ind.keep
      ind.rem <- ind.rem[!pass.thr]
    }

    colnames(scores.all) <- thr.list
  }
  rownames(scores.all) <- ind.test

  scores.all
}

################################################################################

#' Thresholding and correction
#'
#' P-value thresholding and correction of summary statistics for winner's curse.
#'
#' @param beta Vector of effect sizes.
#' @param beta_se Vector of standard errors for `beta`.
#'   Either `beta_se` or `lpS` must be provided.
#' @param lpS Vector of -log10(p-value) associated with `beta`.
#'   Either `beta_se` or `lpS` must be provided.
#' @param thr_lpS Threshold on `lpS` (-log10(p-value) at which variants are
#'   excluded if they  not significant enough.
#'
#' @return `beta` after p-value thresholding and shrinkage.
#' @export
#'
#' @references
#' Zhong, H., & Prentice, R. L. (2008). Bias-reduced estimators and confidence
#' intervals for odds ratios in genome-wide association studies.
#' Biostatistics, 9(4), 621-634.
#'
#' @examples
#' beta <- rnorm(1000)
#' beta_se <- runif(1000, min = 0.3, max = 0.5)
#' new_beta <- snp_thr_correct(beta, beta_se = beta_se, thr_lpS = 1)
#' plot(beta / beta_se, new_beta / beta_se, pch = 20); abline(0, 1, col = "red")
#' plot(beta, new_beta, pch = 20); abline(0, 1, col = "red")
#'
#' # Can provide -log10(p-values) instead of standard errors
#' lpval <- -log10(pchisq((beta / beta_se)^2, df = 1, lower.tail = FALSE))
#' new_beta2 <- snp_thr_correct(beta, lpS = lpval, thr_lpS = 1)
#' all.equal(new_beta2, new_beta)
#'
snp_thr_correct <- function(beta, beta_se, lpS, thr_lpS) {

  if (thr_lpS < 0) stop2("'thr_lpS' must be positive (or 0).")
  if (thr_lpS == 0) return(beta)

  if (!missing(beta_se)) {
    bigassertr::assert_lengths(beta, beta_se)
    z <- abs(beta / beta_se)
  } else if (!missing(lpS)) {
    bigassertr::assert_lengths(beta, lpS)
    z <- sqrt(stats::qchisq(-lpS / log10(exp(1)), log.p = TRUE,
                            df = 1, lower.tail = FALSE))
  } else {
    stop2("'beta_se' and 'lpS' cannot be both missing.")
  }

  thr_Z <- sqrt(stats::qchisq(10^-thr_lpS, df = 1, lower.tail = FALSE))
  Z <- seq(0, 10 * max(z), length.out = 1e6)
  Z2 <- Z + (stats::dnorm(Z - thr_Z) - stats::dnorm(-Z - thr_Z)) /
    (stats::pnorm(Z - thr_Z) + stats::pnorm(-Z - thr_Z))
  knn <- bigutilsr::knn_parallel(Z2, as.matrix(z), k = 1, ncores = 1)
  new_z <- Z[drop(knn$nn.idx)]

  ifelse(z >= thr_Z, beta * pmin(new_z / z, 1), 0)
}

################################################################################
