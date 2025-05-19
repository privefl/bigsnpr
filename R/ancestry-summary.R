################################################################################

#' Estimation of ancestry proportions
#'
#' Estimation of ancestry proportions. Make sure to match summary statistics
#' using [snp_match()] (and to reverse frequencies correspondingly).
#'
#' @param freq Vector of frequencies from which to estimate ancestry proportions.
#' @param info_freq_ref A data frame (or matrix) with the set of frequencies to
#'   be used as reference (one population per column).
#' @param projection Matrix of "loadings" for each variant/PC to be used to
#'   project allele frequencies.
#' @param correction Coefficients to correct for shrinkage when projecting.
#' @param min_cor Minimum correlation between observed and predicted frequencies.
#'   Default is 0.4. When correlation is lower, an error is returned.
#'   For individual genotypes, this should be larger than 0.6.
#'   For allele frequencies, this should be larger than 0.9.
#' @param sum_to_one Whether to force ancestry coefficients to sum to 1?
#'   Default is `TRUE` (otherwise, the sum can be lower than 1).
#'
#' @return Vector of coefficients representing the ancestry proportions.
#'   Also (as attributes) `cor_each`, the correlation between input
#'   frequencies and each reference frequencies, and `cor_pred`, the correlation
#'   between input and predicted frequencies.
#' @export
#'
#' @importFrom stats cor
#'
#' @example examples/example-ancestry-summary.R
#'
snp_ancestry_summary <- function(freq, info_freq_ref, projection, correction,
                                 min_cor = 0.4, sum_to_one = TRUE) {

  assert_package("quadprog")
  assert_nona(freq)
  assert_nona(info_freq_ref)
  assert_nona(projection)
  assert_lengths(freq, rows_along(info_freq_ref), rows_along(projection))
  assert_lengths(correction, cols_along(projection))

  X0 <- as.matrix(info_freq_ref)
  if (mean(cor(X0, freq)) < -0.2)
    stop2("Frequencies seem all reversed; switch reference allele?")

  # project allele frequencies onto the PCA space
  projection <- as.matrix(projection)
  X <- crossprod(projection, X0)
  y <- crossprod(projection, freq) * correction

  cp_X_pd <- Matrix::nearPD(crossprod(X), base.matrix = TRUE)
  if (!isTRUE(cp_X_pd$converged))
    stop2("Could not find nearest positive definite matrix.")

  # solve QP problem using https://stats.stackexchange.com/a/21566/135793
  res <- quadprog::solve.QP(
    Dmat = cp_X_pd$mat,
    dvec = crossprod(y, X),
    Amat = cbind(-1, diag(ncol(X))),
    bvec = c(-1, rep(0, ncol(X))),
    meq  = `if`(sum_to_one, 1, 0)
  )

  cor_pred <- drop(cor(drop(X0 %*% res$solution), freq))
  if (cor_pred < min_cor)
    stop2("Correlation between frequencies is too low: %.3f; %s",
          cor_pred, "check matching between variants.")
  if (cor_pred < 0.99)
    warning2("The solution does not perfectly match the frequencies.")

  structure(round(res$solution, 7),
            names = colnames(info_freq_ref),
            cor_each = drop(cor(X0, freq)),
            cor_pred = cor_pred)
}

################################################################################
