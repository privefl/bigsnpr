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
#'
#' @return vector of coefficients representing the ancestry proportions.
#' @export
#'
#' @importFrom stats cor
#'
#' @example examples/example-ancestry-summary.R
#'
snp_ancestry_summary <- function(freq, info_freq_ref, projection, correction) {

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
    Amat = cbind(1, diag(ncol(X))),
    bvec = c(1, rep(0, ncol(X))),
    meq  = 1
  )

  cor_pred <- drop(cor(drop(X0 %*% res$solution), freq))
  if (cor_pred < 0.99)
    warning2("The solution does not perfectly match the frequencies.")

  setNames(round(res$solution, 7), colnames(info_freq_ref))
}

################################################################################
