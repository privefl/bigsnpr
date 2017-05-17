################################################################################

#' PRS
#'
#' Polygenic Risk Scores with possible clumping and thresholding.
#'
#' @inheritParams bigsnpr-package
#' @param betas Numeric vector of weights associated with each SNP.
#' You may want to see [big_univLinReg] or [big_univLogReg].
#' @param ind.test The individuals on whom to project the scores.
#' @param ind.keep Column (SNP) indices to use (if using clumping, the
#' output of [snp_clumping]). Default doesn't clump.
#' @param lpS Numeric vector of `-log10(p.value)` associated with `betas`.
#' Default doesn't use thresholding.
#' @param thr.list Threshold vector on `lpS` at which SNPs are excluded if
#' they are not significant enough. Default doesn't use thresholding.
#'
#' @return A matrix of scores, where rows correspond to `ind.test` and
#' columns correspond to `thr.list`.
#' @export
#'
#' @example examples/example-PRS.R
snp_PRS <- function(G, betas, ind.test, ind.keep = cols_along(G),
                    lpS = NULL, thr.list = 0) {

  check_args()

  # thresholding and projecting
  if (is.null(lpS) || isTRUE(all.equal(thr.list, 0))) {
    message("'lpS' or 'thr.list' was not specified. Thresholding disabled.")
    scores.all <- as.matrix(big_prodVec(G, betas[ind.keep], ind.test, ind.keep))
  } else {
    n.thr <- length(thr.list)

    scores.all <- foreach(ic = 1:n.thr, .combine = 'cbind') %do% {

      ind.col <- intersect(ind.keep, which(lpS > thr.list[ic]))
      if (length(ind.col)) {
        big_prodVec(G, betas[ind.col], ind.test, ind.col)
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
