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

imputeChr <- function(Gna, infos.imp, ind.chr, alpha, size, p.train, n.cor, seed,
                      thr.imp) {

  # Do something only if there is something to do
  need_imp <- is.na(infos.imp[1, ind.chr])
  if (any(need_imp)) {

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

    # Impute (using mode) some SNPs with few missing values or low variation
    if (!is.null(thr.imp)) {
      counts <- big_counts(Gna, ind.col = ind.chr)
      which_max <- max.col(t(counts[1:3, ]))
      maxs <- counts[cbind(which_max, cols_along(counts))]
      nbNA <- counts[4, ]
      need_imp <- need_imp & (nbNA > 0)
      infos.imp[1, ind.chr[nbNA == 0]] <- 0
      pNA <- nbNA / n
      pDiff <- (n + 1 - maxs) / (n + 1 - nbNA)
      can_imp_by_mode <- ((pNA * pDiff) < thr.imp)
      imp_by_mode <- (need_imp & can_imp_by_mode)
      big_apply(Gna, function(X, ind) {
        ind2 <- ind.chr[ind]
        indNA <- which(is.na(X[, ind2, drop = FALSE]), arr.ind = TRUE)
        X[cbind(indNA[, 1], ind2[indNA[, 2]])] <- (which_max + 3L)[ind[indNA[, 2]]]
        NULL
      }, ind = which(imp_by_mode))
      infos.imp[, ind.chr[imp_by_mode]] <- rbind(pNA, pDiff)[, imp_by_mode]
      need_imp <- (need_imp & !can_imp_by_mode)
    }

    # imputation
    for (i in which(need_imp)) {

      snp <- ind.chr[i]
      X.label <- X[, snp]
      l <- length(indNA <- which(is.na(X.label)))
      if (l > 0) {
        indNoNA <- setdiff(seq_len(n), indNA)
        ind.train <- sort(sample(indNoNA, size = p.train * length(indNoNA)))
        ind.val <- setdiff(indNoNA, ind.train)

        ind.col <- ind.chr[which(corr[, i] != 0)]
        if (length(ind.col) < 5L)
          ind.col <- intersect(setdiff(-size:size + snp, snp), ind.chr)

        bst <- xgboost(
          data = X2[ind.train, ind.col, drop = FALSE],
          label = X.label[ind.train],
          objective = "binary:logistic",
          base_score = mean(X.label[ind.train]) / 2,
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
#'   Default uses maximum 10K individuals for estimating correlations.
#' @param seed An integer, for reproducibility. Default doesn't use seeds.
#' @param thr.imp Default is `NULL` and do nothing. Otherwise, you can use a
#'   proportion (e.g. `0.01`), which is the maximum estimated proportion of
#'   wrongly imputed individuals if imputing using the mode by SNP.
#'
#' @return An [FBM] with
#' - the proportion of missing values by SNP (first row),
#' - the estimated proportion of imputation errors by SNP (second row).
#' @export
#'
#' @import Matrix xgboost
#'
#' @example examples/example-impute.R
#'
snp_fastImpute <- function(Gna, infos.chr,
                           alpha = 1e-4,
                           size = 200,
                           p.train = 0.8,
                           n.cor = min(nrow(Gna), 10e3),
                           seed = NA,
                           ncores = 1,
                           thr.imp = NULL) {

  check_args(infos.chr = "assert_lengths(infos.chr, cols_along(Gna))")

  infos.imp <- FBM_infos(Gna)

  if (!is.na(seed)) seed <- seq_len(max(infos.chr)) + seed
  args <- as.list(environment())

  do.call(what = snp_split, args = c(args, FUN = imputeChr, combine = 'c'))

  infos.imp
}

################################################################################
