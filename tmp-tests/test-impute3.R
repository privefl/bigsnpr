################################################################################

CODE_IMPUTE_LABEL <- c(0, 0.5, 1, rep(NA, 253))
CODE_IMPUTE_PRED  <- c(0, 1, 2, NA, 0, 1, 2, rep(NA, 249))

################################################################################

imputeChr <- function(Gna, ind.chr, alpha, size, p.train, seed) {

  # reproducibility
  if (!any(is.na(seed))) set.seed(seed[attr(ind.chr, "chr")])

  # init
  n <- nrow(Gna)
  m.chr <- length(ind.chr)

  loss_fun <- function(x, y, t1 = 0.25, t2 = 0.75, lambda = 0) {
    mean(((x > t1) + (x > t2) - y)^2) +
      lambda * ((t1 - 0.25)^2 + (t2 - 0.75)^2)
  }

  # correlation between SNPs
  corr <- snp_cor(
    Gna = Gna,
    ind.row = 1:n,
    ind.col = ind.chr,
    size = size,
    alpha = alpha,
    fill.diag = FALSE
  )

  # imputation
  nbNA <- integer(m.chr)
  error <- rep(NA_real_, m.chr)
  num_pred <- rep(NA_integer_, m.chr)
  for (i in 1:m.chr) {
    cat(i)
    X.label <- Gna[, ind.chr[i]]
    nbNA[i] <- l <- length(indNA <- which(is.na(X.label)))
    if (l > 0) {
      indNoNA <- setdiff(1:n, indNA)
      ind.train <- sort(sample(indNoNA, p.train * length(indNoNA)))
      ind.val <- setdiff(indNoNA, ind.train)

      # ind.col <- ind.chr[which(corr[, i] != 0)]
      ind.col <- which(corr[, i] != 0)
      num_pred[i] <- length(ind.col)
      if (length(ind.col) < 5)
        ind.col <- setdiff(intersect(1:m.chr, -size:size + i), i)
      ind.col <- ind.chr[ind.col]

      # xgboost model
      bst <- xgboost(data = Gna[ind.train, ind.col],
                     label = X.label[ind.train] / 2,
                     objective = "binary:logistic",
                     base_score = mean(X.label[ind.train]) / 2,
                     nrounds = 10,
                     params = list(max_depth = 4, gamma = 1, alpha = 1),
                     nthread = 1, verbose = 0, save_period = NULL)
      # learn thresholds on training set
      pred.train <- stats::predict(bst, Gna[ind.train, ind.col])
      lambda <- 8 * loss_fun(pred.train, X.label[ind.train])
      opt.min <- stats::optim(par = c(0.25, 0.75), fn = function(t) {
        loss_fun(pred.train, X.label[ind.train], t[[1]], t[[2]], lambda)
      })
      thrs <- `if`(opt.min$convergence == 0, opt.min$par, c(0.25, 0.75))

      # error of validation
      pred.val <- stats::predict(bst, Gna[ind.val, ind.col, drop = FALSE])
      pred.val <- rowSums(outer(pred.val, thrs, '>'))
      error[i] <- mean(pred.val != X.label[ind.val])
      # impute
      pred <- stats::predict(bst, Gna[indNA, ind.col, drop = FALSE])
      pred <- rowSums(outer(pred, thrs, '>'))
      Gna[indNA, ind.chr[i]] <- as.raw(pred + 4)
    }
  }

  data.frame(pNA = nbNA / n, pError = error, num_pred = num_pred)
}

################################################################################

#' Fast imputation
#'
#' Fast imputation algorithm based on local XGBoost models. **This algorithm
#' has not been extensively compared with other imputation methods yet.**
#'
#' @inheritParams bigsnpr-package
#' @param alpha Type-I error for testing correlations. Default is `0.02`.
#' @param size Number of neighbor SNPs to be possibly included in the model
#' imputing this particular SNP. Default is `500`.
#' @param p.train Proportion of non missing genotypes that are used for training
#' the imputation model while the rest is used to assess the accuracy of
#' this imputation model. Default is `0.8`.
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
                           alpha = 0.02,
                           size = 500,
                           p.train = 0.8,
                           seed = NA,
                           ncores = 1) {

  check_args()

  Gna$code256 <- CODE_IMPUTE_PRED

  if (!is.na(seed)) seed <- seq_len(max(infos.chr)) + seed
  args <- as.list(environment())

  do.call(what = snp_split, args = c(args, FUN = imputeChr, combine = 'rbind'))
}

################################################################################

#' imputeChr2 <- function(Gna, ind.chr, size, p.train, seed) {
#'
#'   # reproducibility
#'   if (!any(is.na(seed))) set.seed(seed[attr(ind.chr, "chr")])
#'
#'   # init
#'   n <- nrow(Gna)
#'   m.chr <- length(ind.chr)
#'
#'   # imputation
#'   nbNA <- integer(m.chr)
#'   error <- rep(NA_real_, m.chr)
#'   for (i in 1:m.chr) {
#'     cat(i)
#'     X.label <- Gna[, ind.chr[i]]
#'     nbNA[i] <- l <- length(indNA <- which(is.na(X.label)))
#'     if (l > 0) {
#'       indNoNA <- setdiff(1:n, indNA)
#'       ind.train <- sort(sample(indNoNA, p.train * length(indNoNA)))
#'       ind.val <- setdiff(indNoNA, ind.train)
#'
#'       ind.col <- -size:size + i
#'       ind.col[ind.col < 1 | ind.col > m.chr | ind.col == i] <- 0L
#'       X.data <- Gna[, ind.col]
#'
#'       bst <- xgboost(
#'         data = X.data[ind.train, ],
#'         label = X.label[ind.train],
#'         objective = "multi:softmax",
#'         base_score = mean(X.label[ind.train]),
#'         nrounds = 10,
#'         params = list(max_depth = 4, num_class = 3, gamma = 1, alpha = 1),
#'         nthread = 1,
#'         verbose = 0,
#'         save_period = NULL
#'       )
#'
#'       # error of validation
#'       pred.val <- stats::predict(bst, X.data[ind.val, ])
#'       error[i] <- mean(pred.val != X.label[ind.val])
#'       # impute
#'       pred <- stats::predict(bst, X.data[indNA, ])
#'       Gna[indNA, ind.chr[i]] <- as.raw(pred + 4L)
#'     }
#'   }
#'
#'   data.frame(pNA = nbNA / n, pError = error)
#' }
#'
#' ################################################################################
#'
#' #' Fast imputation
#' #'
#' #' Fast imputation algorithm based on local XGBoost models. **This algorithm
#' #' has not been extensively compared with other imputation methods yet.**
#' #'
#' #' @inheritParams bigsnpr-package
#' #' @param size Number of neighbor SNPs to be possibly included in the model
#' #' imputing this particular SNP. Default is `100`.
#' #' @param p.train Proportion of non missing genotypes that are used for training
#' #' the imputation model while the rest is used to assess the accuracy of
#' #' this imputation model. Default is `0.8`.
#' #' @param seed An integer, for reproducibility. Default doesn't use seeds.
#' #'
#' #' @return A `data.frame` with
#' #' - the proportion of missing values by SNP,
#' #' - the estimated proportion of imputation errors by SNP.
#' #' @export
#' #'
#' #' @import xgboost
#' #'
#' snp_fastImpute2 <- function(Gna, infos.chr,
#'                            size = 100,
#'                            p.train = 0.8,
#'                            seed = NA,
#'                            ncores = 1) {
#'
#'   check_args()
#'
#'   Gna$code256 <- CODE_IMPUTE_PRED
#'
#'   if (!is.na(seed)) seed <- seq_len(max(infos.chr)) + seed
#'   args <- as.list(environment())
#'
#'   do.call(what = snp_split, args = c(args, FUN = imputeChr2, combine = 'rbind'))
#' }
#'
#' ################################################################################
#'
#' imputeChr3 <- function(Gna, ind.chr, alpha, size, p.train, seed) {
#'
#'   # reproducibility
#'   if (!any(is.na(seed))) set.seed(seed[attr(ind.chr, "chr")])
#'
#'   # init
#'   X <- Gna$copy(code = CODE_IMPUTE_PRED)
#'   n <- nrow(X)
#'   m.chr <- length(ind.chr)
#'
#'   # correlation between SNPs
#'   corr <- snp_cor(
#'     Gna = Gna,
#'     ind.row = 1:n,
#'     ind.col = ind.chr,
#'     size = size,
#'     alpha = alpha,
#'     fill.diag = FALSE
#'   )
#'   print(corr)
#'
#'   # imputation
#'   nbNA <- integer(m.chr)
#'   error <- rep(NA_real_, m.chr)
#'   for (i in 1:m.chr) {
#'     cat(i)
#'     X.label <- Gna[, ind.chr[i]]
#'     nbNA[i] <- l <- length(indNA <- which(is.na(X.label)))
#'     if (l > 0) {
#'       indNoNA <- setdiff(1:n, indNA)
#'       ind.train <- sort(sample(indNoNA, p.train * length(indNoNA)))
#'       ind.val <- setdiff(indNoNA, ind.train)
#'
#'       ind.col <- which(corr[, i] != 0)
#'       if (length(ind.col) < 5)
#'         ind.col <- setdiff(intersect(1:m.chr, -size:size + i), i)
#'       X.data <- X[, ind.chr[ind.col]]
#'
#'       bst <- xgboost(
#'         data = X.data[ind.train, ],
#'         label = X.label[ind.train],
#'         objective = "multi:softmax",
#'         base_score = mean(X.label[ind.train]),
#'         nrounds = 10,
#'         params = list(max_depth = 4, num_class = 3, gamma = 1, alpha = 1),
#'         nthread = 1,
#'         verbose = 0,
#'         save_period = NULL
#'       )
#'
#'       # error of validation
#'       pred.val <- stats::predict(bst, X.data[ind.val, ])
#'       error[i] <- mean(pred.val != X.label[ind.val])
#'       # impute
#'       pred <- stats::predict(bst, X.data[indNA, ])
#'       Gna[indNA, ind.chr[i]] <- as.raw(pred + 4L)
#'     }
#'   }
#'
#'   data.frame(pNA = nbNA / n, pError = error)
#' }
#'
#' ################################################################################
#'
#' #' Fast imputation
#' #'
#' #' Fast imputation algorithm based on local XGBoost models. **This algorithm
#' #' has not been extensively compared with other imputation methods yet.**
#' #'
#' #' @inheritParams bigsnpr-package
#' #' @param alpha Type-I error for testing correlations. Default is `0.02`.
#' #' @param size Number of neighbor SNPs to be possibly included in the model
#' #' imputing this particular SNP. Default is `500`.
#' #' @param p.train Proportion of non missing genotypes that are used for training
#' #' the imputation model while the rest is used to assess the accuracy of
#' #' this imputation model. Default is `0.8`.
#' #' @param seed An integer, for reproducibility. Default doesn't use seeds.
#' #'
#' #' @return A `data.frame` with
#' #' - the proportion of missing values by SNP,
#' #' - the estimated proportion of imputation errors by SNP.
#' #' @export
#' #'
#' #' @import Matrix xgboost
#' #'
#' snp_fastImpute3 <- function(Gna, infos.chr,
#'                            alpha = 0.02,
#'                            size = 500,
#'                            p.train = 0.8,
#'                            seed = NA,
#'                            ncores = 1) {
#'
#'   check_args()
#'
#'   if (!is.na(seed)) seed <- seq_len(max(infos.chr)) + seed
#'   args <- as.list(environment())
#'
#'   do.call(what = snp_split, args = c(args, FUN = imputeChr3, combine = 'rbind'))
#' }
#'
#' ################################################################################
