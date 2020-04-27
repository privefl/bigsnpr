################################################################################

CODE_IMPUTE_LABEL <- c(0, 0.5, 1, rep(NA, 253))

#' CODE_IMPUTE_PRED: code genotype calls and missing values (4),
#' and imputed calls (3).
#' @rdname snp_codes
#' @export
CODE_IMPUTE_PRED  <- c(0, 1, 2, NA, 0, 1, 2, rep(NA, 249))

################################################################################

FBM_infos <- function(Gna) {

  base <- sub_bk(Gna$backingfile, "-infos-impute", stop_if_not_ext = FALSE)
  rds <- paste0(base, ".rds")

  if (file.exists(rds)) {
    # cat("Attach\n")
    readRDS(rds)
  } else {
    # cat("Create\n")
    FBM(2, ncol(Gna), backingfile = base, init = NA_real_)$save()
  }
}

################################################################################

imputeChr <- function(X, X2, infos.imp, ind.chr, alpha, size,
                      p.train, n.cor, seed, ncores) {

  old <- .Random.seed
  on.exit(.Random.seed <<- old, add = TRUE)

  # Do something only if there is still something to do
  if (any(is.na(infos.imp[1, ind.chr]))) {

    n <- nrow(X)

    if (!is.na(seed)) set.seed(seed)

    # correlation between variants
    corr <- snp_cor(
      Gna       = X,
      ind.row   = sort(sample(n, size = n.cor)),
      ind.col   = ind.chr,
      size      = size,
      alpha     = alpha,
      fill.diag = FALSE,
      ncores    = ncores
    )

    # imputation
    for (i in seq_along(ind.chr)) {

      if (!is.na(seed)) set.seed(seed + i)

      snp <- ind.chr[i]
      # Do something only if it wasn't done before
      if (is.na(infos.imp[1, snp])) {

        X.label <- X[, snp]
        nbna <- length(indNA <- which(is.na(X.label)))
        if (nbna > 0) {
          indNoNA <- setdiff(seq_len(n), indNA)
          ind.train <- sort(sample(indNoNA, size = p.train * length(indNoNA)))
          ind.val <- setdiff(indNoNA, ind.train)

          ind.col <- ind.chr[which(corr[, i] != 0)]
          if (length(ind.col) < 5L)
            ind.col <- intersect(setdiff(-size:size + snp, snp), ind.chr)

          data.train <- xgboost::xgb.DMatrix(
            label = X.label[ind.train],
            data  = X2[ind.train, ind.col, drop = FALSE])

          bst.params <- list(
            objective  = "binary:logistic",
            max_depth  = 4,
            base_score = min(max(1e-7, mean(X.label[ind.train])), 1 - 1e-7),
            verbose    = 0,
            nthread    = ncores
          )

          bst <- xgboost::xgb.train(
            data    = data.train,
            params  = bst.params,
            nrounds = 10
          )

          # error of validation
          pred2 <- stats::predict(bst, X2[ind.val, ind.col, drop = FALSE])
          infos.imp[2, snp] <- mean(round(2 * pred2) != (2 * X.label[ind.val]))
          # imputation
          pred <- stats::predict(bst, X2[indNA, ind.col, drop = FALSE])
          X2[indNA, snp] <- as.raw(round(2 * pred) + 4)
        }

        # this variant is done
        infos.imp[1, snp] <- nbna / n
      }
    }

  }

  invisible()
}

################################################################################

#' Fast imputation
#'
#' Fast imputation algorithm based on local XGBoost models.
#'
#' @inheritParams bigsnpr-package
#' @param alpha Type-I error for testing correlations. Default is `1e-4`.
#' @param size Number of neighbor SNPs to be possibly included in the model
#'   imputing this particular SNP. Default is `200`.
#' @param p.train Proportion of non missing genotypes that are used for training
#'   the imputation model while the rest is used to assess the accuracy of
#'   this imputation model. Default is `0.8`.
#' @param n.cor Number of rows that are used to estimate correlations.
#'   Default uses them all.
#' @param seed An integer, for reproducibility. Default doesn't use seeds.
#'
#' @return An [FBM] with
#' - the proportion of missing values by SNP (first row),
#' - the estimated proportion of imputation errors by SNP (second row).
#' @export
#'
#' @import Matrix
#'
#' @seealso [snp_fastImputeSimple()]
#'
#' @example examples/example-impute.R
#'
snp_fastImpute <- function(Gna, infos.chr,
                           alpha = 1e-4,
                           size = 200,
                           p.train = 0.8,
                           n.cor = nrow(Gna),
                           seed = NA,
                           ncores = 1) {

  check_args(infos.chr = "assert_lengths(infos.chr, cols_along(Gna))")
  assert_package("xgboost")

  X  <- Gna$copy(code = CODE_IMPUTE_LABEL)
  X2 <- Gna$copy(code = CODE_IMPUTE_PRED)

  infos.imp <- FBM_infos(Gna)

  ind.chrs <- split(seq_along(infos.chr), infos.chr)
  for (ind in ind.chrs) {
    imputeChr(X, X2, infos.imp, ind, alpha, size, p.train, n.cor, seed, ncores)
  }

  infos.imp
}

################################################################################

#' Fast imputation
#'
#' Fast imputation via mode, mean, sampling according to allele frequencies, or 0.
#'
#' @inheritParams bigsnpr-package
#' @param method Either `"random"` (sampling according to allele frequencies),
#'   `"mean0"` (rounded mean), `"mean2"` (rounded mean to 2 decimal places),
#'   `"mode"` (most frequent call).
#'
#' @return A new `FBM.code256` object (same file, but different code).
#' @export
#'
#' @seealso [snp_fastImpute()]
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
#' G$copy(code = c(0, 1, 2, rep(0, 253)))[, 2]  # "imputation" by 0
#'
snp_fastImputeSimple <- function(
  Gna, method = c("mode", "mean0", "mean2", "random"), ncores = 1) {

  check_args()
  stopifnot(identical(Gna$code256, CODE_012))

  if (identical(method, "zero")) {
    warning2("Using 'method = \"zero\"' is deprecated. Using $copy() instead..")
    Gna$copy(code = c(0, 1, 2, 0, rep(NA, 252)))
  } else {
    method <- match(match.arg(method), c("mode", "mean0", "mean2", "random"))
    impute(Gna, method, ncores)
    Gna$copy(code = `if`(method == 3, CODE_DOSAGE, CODE_IMPUTE_PRED))
  }
}
