bigsnp <- snp_attachExtdata("example-missing.bed")
G <- bigsnp$genotypes
G[, 2]
ind <- 2:15
G_sub <- G[, ind, drop = FALSE]

ind_NA <- which(is.na(G_sub), arr.ind = TRUE)
col_NA_rel <- ind_NA[, 2]
col_NA_abs <- ind[col_NA_rel]
means <- colMeans(G_sub, na.rm = TRUE)

## sampling according to frequencies
G[cbind(ind_NA[, 1], col_NA_abs)] <-
  rbinom(length(col_NA_abs), size = 2, prob = (means / 2)[col_NA_rel]) + 4L
G[, 2]
G$copy(code = CODE_IMPUTE_PRED)[, 2]

## rounded mean (0 decimal places)
G[cbind(ind_NA[, 1], col_NA_abs)] <-
  as.raw(round(means) + 4)[col_NA_rel]
G[, 2]
G$copy(code = CODE_IMPUTE_PRED)[, 2]
mean(G[, 2], na.rm = TRUE)

## rounded mean (2 decimal places)
G[cbind(ind_NA[, 1], col_NA_abs)] <-
  as.raw(round(100 * means) + 7)[col_NA_rel]
G[, 2]
G$copy(code = CODE_IMPUTE_PRED)[, 2]
G$copy(code = CODE_DOSAGE)[, 2]
mean(G[, 2], na.rm = TRUE)

## mode
counts <- big_counts(G, ind.col = ind)[1:3, ]
G[cbind(ind_NA[, 1], col_NA_abs)] <-
  as.raw(apply(counts, 2, which.max) + 3L)[col_NA_rel]
G[, 2]
G$copy(code = CODE_IMPUTE_PRED)[, 2]
G$copy(code = CODE_DOSAGE)[, 2]
mean(G[, 2], na.rm = TRUE)

x <- "lol"
tmp <- switch (x,
        blabla = 2:7,
        lol = 1:100
)
tmp



#' Fast imputation
#'
#' Fast imputation via mode, mean or sampling according to allele frequencies.
#'
#' @inheritParams bigsnpr-package
#' @param method Either `"random"` (sampling according to allele frequencies),
#'   `"mean0"` (rounded mean), `"mean2"` (rounded mean to 2 decimal places),
#'   `"mode"` (most frequent call).
#'
#' @return A new `FBM.code256` object (same file, but different code).
#' @export
#'
#' @examples
#' bigsnp <- snp_attachExtdata("example-missing.bed")
#' G <- bigsnp$genotypes
#' G[, 2]  # some missing values
#' G2 <- snp_fastImputeSimple(G)
#' G2[, 2]  # no missing values anymore
#' G[, 2]  # imputed, but still returning missing values
#' G$copy(code = CODE_IMPUTE_PRED)[, 2]  # need to decode imputed values
#'
snp_fastImputeSimple <- function(
  Gna, method = c("random", "mean0", "mean2", "mode"), ncores = 1) {

  check_args()

  stopifnot(identical(Gna$code256, CODE_012))

  method <- match.arg(method)
  CODE <- `if`(method == "mean2", CODE_DOSAGE, CODE_IMPUTE_PRED)

  big_apply(Gna, function(X, ind, method) {

    G_sub <- X[, ind, drop = FALSE]

    ind_NA <- which(is.na(G_sub), arr.ind = TRUE)
    col_NA_rel <- ind_NA[, 2]
    col_NA_abs <- ind[col_NA_rel]

    Gna[cbind(ind_NA[, 1], col_NA_abs)] <- switch(
      method,
      mode = {
        counts <- big_counts(Gna, ind.col = ind)[1:3, ]
        as.raw(apply(counts, 2, which.max) + 3L)[col_NA_rel]
      },
      random = {
        af <- colMeans(G_sub, na.rm = TRUE) / 2
        rbinom(length(col_NA_abs), size = 2, prob = af[col_NA_rel]) + 4L
      },
      mean0 = {
        means <- colMeans(G_sub, na.rm = TRUE)
        as.raw(round(means) + 4)[col_NA_rel]
      },
      mean2 = {
        means <- colMeans(G_sub, na.rm = TRUE)
        as.raw(round(100 * means) + 7)[col_NA_rel]
      }
    )

    NULL
  }, a.combine = 'c', ncores = ncores, method = method)

  Gna$copy(code = CODE)
}
