################################################################################

CODE_IMPUTE_LABEL <- c(0, 0.5, 1, rep(NA, 253))
CODE_IMPUTE_PRED  <- c(0, 1, 2, NA, 0, 1, 2, rep(NA, 249))

################################################################################

imputeChr <- function(Gna, ind.chr, alpha, size, p.train, n.cor, seed) {

  # reproducibility
  if (!any(is.na(seed))) set.seed(seed[attr(ind.chr, "chr")])

  # init
  X <- Gna$copy(code = CODE_IMPUTE_LABEL)
  n <- nrow(X)
  m.chr <- length(ind.chr)
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
  nbNA <- integer(m.chr)
  error <- rep(NA_real_, m.chr)
  for (i in 1:m.chr) {
    snp <- ind.chr[i]
    X.label <- X[, snp]
    nbNA[i] <- l <- length(indNA <- which(is.na(X.label)))
    if (l > 0) {
      indNoNA <- setdiff(1:n, indNA)
      ind.train <- sort(sample(indNoNA, p.train * length(indNoNA)))
      ind.val <- setdiff(indNoNA, ind.train)

      ind.col <- ind.chr[which(corr[, i] != 0)]
      if (length(ind.col) < 5L)
        ind.col <- intersect(setdiff(-size:size + snp, snp), ind.chr)

      bst <- xgboost(data = X2[ind.train, ind.col, drop = FALSE],
                     label = X.label[ind.train],
                     objective = "binary:logistic",
                     base_score = mean(X.label[ind.train]) / 2,
                     nrounds = 10,
                     params = list(max_depth = 4),
                     nthread = 1,
                     verbose = 0,
                     save_period = NULL)

      # error of validation
      pred2 <- stats::predict(bst, X2[ind.val, ind.col, drop = FALSE])
      error[i] <- mean(round(2 * pred2) != (2 * X.label[ind.val]))
      # impute
      pred <- stats::predict(bst, X2[indNA, ind.col, drop = FALSE])
      X2[indNA, snp] <- as.raw(round(2 * pred) + 4)
    }
  }

  data.frame(pNA = nbNA / n, pError = error)
}

################################################################################

#' Fast imputation
#'
#' Fast imputation algorithm based on local XGBoost models. **This algorithm
#' has not been extensively compared with other imputation methods yet.**
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
#' @return A `data.frame` with
#' - the proportion of missing values by SNP,
#' - the estimated proportion of imputation errors by SNP.
#' @export
#'
#' @import Matrix xgboost
#'
snp_fastImpute <- function(Gna, infos.chr,
                           alpha = 1e-4,
                           size = 200,
                           p.train = 0.8,
                           n.cor = nrow(Gna),
                           seed = NA,
                           ncores = 1) {

  check_args(infos.chr = "assert_lengths(infos.chr, cols_along(Gna))")

  if (!is.na(seed)) seed <- seq_len(max(infos.chr)) + seed
  args <- as.list(environment())

  do.call(what = snp_split, args = c(args, FUN = imputeChr, combine = 'rbind'))
}

################################################################################
