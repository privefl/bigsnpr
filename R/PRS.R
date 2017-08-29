################################################################################

#' PRS
#'
#' Polygenic Risk Scores with possible clumping and thresholding.
#'
#' @inheritParams bigsnpr-package
#' @param betas.keep Numeric vector of weights associated with each SNP
#'   corresponding to `ind.keep`.
#'   You may want to see [big_univLinReg] or [big_univLogReg].
#' @param ind.test The individuals on whom to project the scores.
#' @param ind.keep Column (SNP) indices to use (if using clumping, the
#'   output of [snp_clumping]). Default doesn't clump.
#' @param lpS.keep Numeric vector of `-log10(p.value)` associated with
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
snp_PRS <- function(G, betas.keep, ind.test, ind.keep = cols_along(G),
                    lpS.keep = NULL, thr.list = 0) {

  check_args()
  assert_lengths(betas.keep, ind.keep)

  # thresholding and projecting
  if (is.null(lpS.keep) || isTRUE(all.equal(thr.list, 0))) {
    message("'lpS.keep' or 'thr.list' was not specified. Thresholding disabled.")
    scores.all <- as.matrix(big_prodVec(G, betas.keep, ind.test, ind.keep))
  } else {
    assert_lengths(lpS.keep, ind.keep)

    scores.all <- foreach(ic = seq_along(thr.list), .combine = 'cbind') %do% {

      pass.thr <- (lpS.keep > thr.list[ic])
      if (any(pass.thr)) {
        big_prodVec(G, betas.keep[pass.thr], ind.test, ind.keep[pass.thr])
      } else {
        rep(0, length(ind.test))
      }
    }
    colnames(scores.all) <- thr.list
  }
  rownames(scores.all) <- ind.test

  scores.all
}

################################################################################
