################################################################################

getMethodFun <- function(method, y.train) {
  if (method == "auto") {
    l <- length(unique(y.train))
    if (l < 2) stop("'y.train' needs at least two different values.")
    `if`(l > 2, bigstatsr::big_univRegLin, bigstatsr::big_univRegLog)
  } else if (method == "lin") {
    bigstatsr::big_univRegLin
  } else if (method == "log") {
    bigstatsr::big_univRegLog
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
#' @param betas Numeric vector of weights associated with each SNP. You may
#' want to see [big_univRegLin][bigstatsr::big_univRegLin] or
#' [big_univRegLog][bigstatsr::big_univRegLog].
#' @param ind.test The individuals on whom to project the scores.
#' @param lpS Numeric vector of `-log10(p.value)` associated with `betas`.
#' Default doesn't use thresholding.
#' @param thr.list Threshold vector on `lpS` at which SNPs are excluded if
#' they are not significant enough. Default doesn't use thresholding.
#' @inheritDotParams snp_clumping -x -exclude
#'
#' @return A matrix of scores, where rows correspond to `ind.test` and
#' columns correspond to `thr.list`.
#' @export
#'
#' @example
snp_PRS <- function(x, betas, ind.test,
                    lpS = NULL, thr.list = 0,
                    ...) {
  check_x(x)
  X <- x$genotypes

  # clumping
  if (is.null(list(...)$S)) {
    message("'S' was not specified. Clumping disabled.")
  } else {
    # exclude some SNPs that will never enter the model
    # in order to accelerate computations
    exclude <- `if`(is.null(lpS), NULL, which(lpS < min(thr.list)))
    ind.keep <- snp_clumping(x, exclude = exclude, ...)
    print(length(ind.keep))
  }

  # thresholding and projecting
  if (is.null(lpS)) {
    message("'lpS' was not specified. Thresholding disabled.")
    scores.all <- as.matrix(
      bigstatsr::big_prodVec(X, betas[ind.keep], ind.test, ind.keep)
    )
  } else {
    n.thr <- length(thr.list)

    scores.all <- foreach(j = 1:n.thr, .combine = 'cbind') %do% {

      ind.col <- intersect(ind.keep, which(lpS > thr.list[j]))
      if (length(ind.col)) {
        bigstatsr::big_prodVec(X, betas[ind.col], ind.test, ind.col)
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

#' Title
#'
#' @param x
#' @param ind.train
#' @param ind.test
#' @param covar.train
#' @param method
#' @param thr.list
#' @param clumping
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
snp_cvPRS <- function(x, ind.train, ind.test, covar.train = NULL,
                      method = "auto",
                      thr.list = NULL,
                      ncores = 1,
                      K = 10,
                      clumping = TRUE, ...) {
  check_x(x)
  if (length(thr.list) < 2)
    stop("Please provide a longer list of threshold values.")
  y.train <- x$fam$affection[ind.train]

  x2 <- x
  x2$genotypes <- describe(x$genotypes)

  FUN <- getMethodFun(method, y.train)

  n <- length(ind.train)
  ind <- sample(rep_len(1:K, n))

  if (is.seq <- (ncores == 1)) {
    registerDoSEQ()
  } else {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
  }
  res <- foreach(ic = 0:K) %dopar% {
    X <- x2$genotypes <- attach.big.matrix(x2$genotypes)

    if (ic == 0) { # main training
      betas <- FUN(X, y.train, ind.train, covar.train)
      S <- `if`(clumping, abs(betas$estim / betas$std.err), NULL)
      snp_PRS(x2, betas$estim, ind.test,
                        lpS = -log10(betas$p.value),
                        thr.list = thr.list,
                        ncores = 1,
                        S = S,
                        ...)
    } else {
      i.test <- which(ind == ic)
      i.train <- setdiff(1:n, i.test)

      betas <- FUN(X, y.train[i.train], ind.train[i.train],
                   covar.train[i.train, ])
      S <- `if`(clumping, abs(betas$estim / betas$std.err), NULL)
      snp_PRS(x2, betas$estim, ind.train[i.test],
                        lpS = -log10(betas$p.value),
                        thr.list = thr.list,
                        ncores = 1,
                        S = S,
                        ...)
    }
  }
}

################################################################################

#' Title
#'
#' @param x
#' @param ind.train
#' @param ind.test
#' @param covar.train
#' @param method
#' @param thr.list
#' @param clumping
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
snp_cvPRS_fast <- function(x, ind.train, ind.test, covar.train = NULL,
                      method = "auto",
                      thr.list = 0,
                      ncores = 1,
                      K = 10, ...) {
  check_x(x)
  if (length(thr.list) < 2)
    stop("Please provide a longer list of threshold values.")
  y.train <- x$fam$affection[ind.train]

  X.desc <- describe(x$genotypes)

  FUN <- getMethodFun(method, y.train)

  n <- length(ind.train)
  ind <- sample(rep_len(1:K, n))

  # GWAS part
  if (is.seq <- (ncores == 1)) {
    registerDoSEQ()
  } else {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
  }
  betas.all <- foreach(ic = 0:K) %dopar% {
    X <- attach.big.matrix(X.desc)

    if (ic == 0) { # main training
      FUN(X, y.train, ind.train, covar.train)
    } else {
      i.train <- which(ind != ic)
      FUN(X, y.train[i.train], ind.train[i.train], covar.train[i.train, ])
    }
  }
  if (!is.seq) parallel::stopCluster(cl)

  # clumping part -> do it only once. Approx but faster.
  betas.main <- betas.all[[1]]
  S <- abs(betas.main$estim / betas.main$std.err)
  ind.keep <- snp_clumping(x, S = S, ncores = ncores, ...)
  print(length(ind.keep))
  n.thr <- length(thr.list)

  # projecting part
  if (is.seq) {
    registerDoSEQ()
  } else {
    cl2 <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl2)
    on.exit(parallel::stopCluster(cl2))
  }
  res <- foreach(ic = 0:K) %dopar% {
    X <- attach.big.matrix(X.desc)

    tmp <- betas.all[[ic + 1]]
    betas <- tmp$estim
    lpS <- -log10(tmp$p.value)
    tmp.ind.test <- `if`(ic == 0, ind.test, ind.train[ind == ic])

    scores.all <- foreach(j = 1:n.thr, .combine = 'cbind') %do% {
      ind.col <- intersect(ind.keep, which(lpS > thr.list[j]))

      if (length(ind.col)) {
        bigstatsr::big_prodVec(X, betas[ind.col], tmp.ind.test, ind.col)
      } else {
        0
      }
    }

    rownames(scores.all) <- tmp.ind.test
    colnames(scores.all) <- thr.list
    scores.all
  }
}

################################################################################
