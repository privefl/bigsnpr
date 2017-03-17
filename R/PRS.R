################################################################################

getMethodFun <- function(method, y.train) {
  if (method == "auto") {
    l <- length(unique(y.train))
    if (l < 2) stop("'y.train' needs at least two different values.")
    `if`(l > 2, big_univRegLin, big_univRegLog)
  } else if (method == "lin") {
    big_univRegLin
  } else if (method == "log") {
    big_univRegLog
  } else {
    stop("Invalid 'method' argument.")
  }
}

################################################################################

#' PRS
#'
#' Polygenic Risk Scores with possible clumping and thresholding.
#'
#' @inheritParams bigsnpr-package
#' @param betas Numeric vector of weights associated with each SNP.
#' You may want to see [big_univRegLin] or [big_univRegLog].
#' @param ind.test The individuals on whom to project the scores.
#' @param ind.keep Column (SNP) indices to use (if using clumping, the
#' output of [snp_clumping]). Default doesn't clump.
#' @param lpS Numeric vector of `-log10(p.value)` associated with `betas`.
#' Default doesn't use thresholding.
#' @param thr.list Threshold vector on `lpS` at which SNPs are excluded if
#' they are not significant enough. Default doesn't use thresholding.
#'
#' @return A matrix of scores, where rows correspond to `ind.test` and
#' columns correspond to `thr.list`.
#' @export
#'
#' @example
snp_PRS <- function(G, betas, ind.test, ind.keep = cols_along(G),
                    lpS = NULL, thr.list = 0) {

  # thresholding and projecting
  if (is.null(lpS)) {
    message("'lpS' was not specified. Thresholding disabled.")
    scores.all <- as.matrix(
      big_prodVec(X, betas[ind.keep], ind.test, ind.keep)
    )
  } else {
    n.thr <- length(thr.list)

    scores.all <- foreach(j = 1:n.thr, .combine = 'cbind') %do% {

      ind.col <- intersect(ind.keep, which(lpS > thr.list[j]))
      if (length(ind.col)) {
        big_prodVec(X, betas[ind.col], ind.test, ind.col)
      } else {
        0
      }
    }
    colnames(scores.all) <- thr.list
  }
  rownames(scores.all) <- ind.test

  scores.all
}

################################################################################

# #' Title
# #'
# #' @param x
# #' @param ind.train
# #' @param ind.test
# #' @param covar.train
# #' @param method
# #' @param thr.list
# #' @param clumping
# #' @param ...
# #'
# #' @return
# #' @export
# #'
# #' @examples
# snp_cvPRS <- function(x, ind.train, ind.test, covar.train = NULL,
#                       method = "auto",
#                       thr.list = NULL,
#                       ncores = 1,
#                       K = 10,
#                       clumping = TRUE, ...) {
#   check_x(x)
#   if (length(thr.list) < 2)
#     stop("Please provide a longer list of threshold values.")
#   y.train <- x$fam$affection[ind.train]
#
#   x2 <- x
#   x2$genotypes <- describe(x$genotypes)
#
#   FUN <- getMethodFun(method, y.train)
#
#   n <- length(ind.train)
#   ind <- sample(rep_len(1:K, n))
#
#   if (is.seq <- (ncores == 1)) {
#     registerDoSEQ()
#   } else {
#     cl <- parallel::makeCluster(ncores)
#     doParallel::registerDoParallel(cl)
#     on.exit(parallel::stopCluster(cl))
#   }
#   res <- foreach(ic = 0:K) %dopar% {
#     X <- x2$genotypes <- attach.big.matrix(x2$genotypes)
#
#     if (ic == 0) { # main training
#       betas <- FUN(X, y.train, ind.train, covar.train)
#       S <- `if`(clumping, abs(betas$estim / betas$std.err), NULL)
#       snp_PRS(x2, betas$estim, ind.test,
#               lpS = -log10(betas$p.value),
#               thr.list = thr.list,
#               ncores = 1,
#               S = S,
#               ...)
#     } else {
#       i.test <- which(ind == ic)
#       i.train <- setdiff(1:n, i.test)
#
#       betas <- FUN(X, y.train[i.train], ind.train[i.train],
#                    covar.train[i.train, ])
#       S <- `if`(clumping, abs(betas$estim / betas$std.err), NULL)
#       snp_PRS(x2, betas$estim, ind.train[i.test],
#               lpS = -log10(betas$p.value),
#               thr.list = thr.list,
#               ncores = 1,
#               S = S,
#               ...)
#     }
#   }
# }
#
# ################################################################################
