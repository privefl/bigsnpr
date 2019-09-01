################################################################################

CODE_IMPUTE_LABEL <- c(0, 0.5, 1, rep(NA, 253))

#' CODE_IMPUTE_PRED: code genotype calls and missing values (4),
#' and imputed calls (3).
#' @rdname snp_codes
#' @export
CODE_IMPUTE_PRED  <- c(0, 1, 2, NA, 0, 1, 2, rep(NA, 249))

################################################################################

FBM_infos <- function(Gna) {

  bk <- Gna$backingfile
  base <- sub("\\.bk$", "-infos-impute", bk)
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

imputeChr <- function(Gna, infos.imp, ind.chr, alpha, size, p.train, n.cor, seed) {

  # Do something only if there is something to do
  if (any(is.na(infos.imp[1, ind.chr]))) {

    # reproducibility
    if (!any(is.na(seed))) set.seed(seed[attr(ind.chr, "chr")])

    # init
    n  <- nrow(Gna)
    X  <- Gna$copy(code = CODE_IMPUTE_LABEL)
    X2 <- Gna$copy(code = CODE_IMPUTE_PRED)

    # correlation between SNPs
    corr <- snp_cor(
      Gna = X,
      ind.row = sort(sample(n, size = n.cor)),
      ind.col = ind.chr,
      size = size,
      alpha = alpha,
      fill.diag = FALSE
    )

    # imputation
    for (i in seq_along(ind.chr)) {

      snp <- ind.chr[i]
      # Do something only if it wasn't done before
      if (is.na(infos.imp[1, snp])) {

        X.label <- X[, snp]
        l <- length(indNA <- which(is.na(X.label)))
        if (l > 0) {
          indNoNA <- setdiff(seq_len(n), indNA)
          ind.train <- sort(sample(indNoNA, size = p.train * length(indNoNA)))
          ind.val <- setdiff(indNoNA, ind.train)

          ind.col <- ind.chr[which(corr[, i] != 0)]
          if (length(ind.col) < 5L)
            ind.col <- intersect(setdiff(-size:size + snp, snp), ind.chr)

          bst <- xgboost::xgboost(
            data = X2[ind.train, ind.col, drop = FALSE],
            label = X.label[ind.train],
            objective = "binary:logistic",
            base_score = min(max(1e-7, mean(X.label[ind.train])), 1 - 1e-7),
            nrounds = 10,
            params = list(max_depth = 4),
            nthread = 1,
            verbose = 0,
            save_period = NULL
          )

          # error of validation
          pred2 <- stats::predict(bst, X2[ind.val, ind.col, drop = FALSE])
          infos.imp[2, snp] <- mean(round(2 * pred2) != (2 * X.label[ind.val]))
          # imputation
          pred <- stats::predict(bst, X2[indNA, ind.col, drop = FALSE])
          X2[indNA, snp] <- as.raw(round(2 * pred) + 4)
        }
        # this SNP is done
        infos.imp[1, snp] <- l / n
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
#' @example examples/example-impute.R
#'
snp_fastImpute <- function(Gna, infos.chr,
                           alpha = 1e-4,
                           size = 200,
                           p.train = 0.8,
                           n.cor = nrow(Gna),
                           seed = NA,
                           ncores = 1) {

  if (!requireNamespace("xgboost", quietly = TRUE))
    stop2("Please install package 'xgboost'.")

  check_args(infos.chr = "assert_lengths(infos.chr, cols_along(Gna))")

  infos.imp <- FBM_infos(Gna)

  if (!is.na(seed)) seed <- seq_len(max(infos.chr)) + seed
  args <- as.list(environment())

  do.call(what = snp_split, args = c(args, FUN = imputeChr, combine = 'c'))

  infos.imp
}

################################################################################

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
  Gna, method = c("mode", "mean0", "mean2", "random"), ncores = 1) {

  check_args()

  stopifnot(identical(Gna$code256, CODE_012))

  method <- match(match.arg(method), c("mode", "mean0", "mean2", "random"))
  big_parallelize(Gna, function(X, ind, method) {
    impute(X, rows_along(X), ind, method)
  }, ncores = ncores, method = method)

  CODE <- `if`(method == 3, CODE_DOSAGE, CODE_IMPUTE_PRED)
  Gna$copy(code = CODE)
}
