################################################################################

#' Determine reference divergence
#'
#' Determine reference divergence while accounting for strand flips.
#'
#' @param ref1 The reference alleles of the first dataset.
#' @param alt1 The alternative alleles of the first dataset.
#' @param ref2 The reference alleles of the second dataset.
#' @param alt2 The alternative alleles of the second dataset.
#'
#' @return A logical vector whether the references alleles are the same.
#'   Missing values can result from missing values in the inputs or from
#'   ambiguous strands.
#' @export
#'
#' @importFrom dplyr %>%
#'
#' @examples
#' same_ref(ref1 = c("A", "C", "T", "G", NA),
#'          alt1 = c("C", "T", "C", "A", "A"),
#'          ref2 = c("A", "C", "A", "A", "C"),
#'          alt2 = c("C", "G", "G", "G", "A"))
same_ref <- function(ref1, alt1, ref2, alt2) {

  ACTG <- c("A", "C", "T", "G")
  REV_ACTG <- stats::setNames(c("T", "G", "A", "C"), ACTG)

  decoder <- expand.grid(list(ACTG, ACTG, ACTG, ACTG)) %>%
    dplyr::mutate(status = dplyr::case_when(
      # BAD: same reference/alternative alleles in a dataset
      (Var1 == Var2) | (Var3 == Var4) ~ NA,
      # GOOD/TRUE: same reference/alternative alleles between datasets
      (Var1 == Var3) & (Var2 == Var4) ~ TRUE,
      # GOOD/FALSE: reverse reference/alternative alleles
      (Var1 == Var4) & (Var2 == Var3) ~ FALSE,
      # GOOD/TRUE: same reference/alternative alleles after strand flip
      (REV_ACTG[Var1] == Var3) & (REV_ACTG[Var2] == Var4) ~ TRUE,
      # GOOD/FALSE: reverse reference/alternative alleles after strand flip
      (REV_ACTG[Var1] == Var4) & (REV_ACTG[Var2] == Var3) ~ FALSE,
      # BAD: the rest
      TRUE ~ NA
    )) %>%
    reshape2::acast(Var1 ~ Var2 ~ Var3 ~ Var4, value.var = "status")

  decoder[cbind(ref1, alt1, ref2, alt2)]
}

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
snp_PRS <- function(G, betas.keep,
                    ind.test = rows_along(G),
                    ind.keep = cols_along(G),
                    same.keep = rep(TRUE, length(ind.keep)),
                    lpS.keep = NULL,
                    thr.list = 0) {

  check_args()
  assert_lengths(betas.keep, ind.keep)
  assert_lengths(same.keep, ind.keep)
  stopifnot(sum(is.na(same.keep)) == 0)
  assert_type(same.keep, "logical")

  # thresholding and projecting
  if (is.null(lpS.keep) || identical(thr.list, 0)) {
    message("'lpS.keep' or 'thr.list' was not specified. Thresholding disabled.")
    scores.all <- as.matrix(prodVecRev(G, betas.keep, same.keep,
                                       ind.test, ind.keep))
  } else {
    assert_lengths(lpS.keep, ind.keep)

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
