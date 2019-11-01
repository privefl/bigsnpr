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
#'   corresponding to `ind.keep`.
#'   You may want to see [big_univLinReg] or [big_univLogReg].
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
