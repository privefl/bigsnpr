### Leads to overfitting

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
