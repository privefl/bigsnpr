################################################################################

#' PRS
#'
#' Polygenic Risk Scores with possible clumping and thresholding.
#'
#' @inheritParams bigsnpr-package
#' @param betas Numeric vector of weights associated with each SNP. You may
#' want to see [big_univRegLin][bigstatsr::big_univRegLin] or
#' [big_univRegLog][bigstatsr::big_univRegLog].
#' @param ind.test The individuals on whom to project the scores.
#' @param lpS Numeric vector of `-log10(p.value)` associated with `betas`.
#' Default doesn't use thresholding.
#' @param thr.list Threshold vector on `lpS` at which SNPs are excluded if
#' they are not significant enough. Default doesn't use thresholding.
#' @inheritDotParams snp_clumping -x -exclude
#'
#' @return A matrix of scores, where rows correspond to `ind.test` and
#' columns correspond to `thr.list`.
#' @export
#'
#' @example
snp_PRS <- function(x, betas, ind.test,
                    lpS = NULL, thr.list = 0,
                    ...) {
  check_x(x)
  X <- x$genotypes

  # clumping
  if (is.null(list(...)$S)) {
    message("'S' was not specified. Clumping disabled.")
  } else {
    # exclude some SNPs that will never enter the model
    # in order to accelerate computations
    exclude <- `if`(is.null(lpS), NULL, which(lpS < min(thr.list)))
    ind.keep <- snp_clumping(x, exclude = exclude, ...)
    print(length(ind.keep))
  }

  # thresholding and projecting
  if (is.null(lpS)) {
    message("'lpS' was not specified. Thresholding disabled.")
    scores.all <- as.matrix(
      bigstatsr::big_prodVec(X, betas[ind.keep], ind.test, ind.keep)
    )
  } else {
    n.thr <- length(thr.list)

    scores.all <- foreach(j = 1:n.thr, .combine = 'cbind') %do% {

      ind.col <- intersect(ind.keep, which(lpS > thr.list[j]))
      if (length(ind.col)) {
        bigstatsr::big_prodVec(X, betas[ind.col], ind.test, ind.col)
      } else {
        0
      }
    }
    colnames(scores.all) <- as.character(thr.list)
  }

  scores.all
}

################################################################################
