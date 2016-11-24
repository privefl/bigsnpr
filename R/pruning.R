#' @title LD pruning for a "bigSNP"
#' @description LD pruning for a \code{bigSNP}.
#' @name Prune
#' @inheritParams bigsnpr-package
#' @param S Numeric vector of summary statistics computed
#' only with `ind.train`.
#' @param pS Numeric vector of p-values associated with S.
#' **`pS` needs to be computed through a decreasing function of `S`**.
#' For example, if `S` follows the standard normal distribution,
#' you should use `abs(S)` instead and compute
#' `pS = 2 * pnorm(abs(S), lower.tail = FALSE)`.
#' @param size Radius of the window's size for the LD evaluations.
#' @param thr.pvalue Threshold on \eqn{-log_{10}(p-value)} to assess
#' which SNPs are kept. Here, it has the purpose to accelerate computations.
#' Default is 1.
#' @param thr.corr Threshold on the correlation between two SNPs.
#' SNPs which are too correlated with another SNP which is more correlated
#' with the disease are pruned.
#' @example examples/example.pruning.R
#' @export
Prune <- function(x,
                  ind.train = seq(nrow(X)),
                  fun.stats,
                  thr.pvalue = 1,
                  size = 2000,
                  thr.corr = 0.2,
                  ncores = 1) {
  #check_x(x, check.y = TRUE)

  # get descriptors
  X <- x$genotypes
  X.desc <- describe(X)

  # export fun.stats and backingpath (force eval of promise)
  PATH <- x$backingpath
  #FUN <- fun.stats


  range.chr <- LimsChr(x)

  if (is.seq <- (ncores == 1)) {
    registerDoSEQ()
  } else {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
  }
  res <- foreach(i = seq_len(nrow(range.chr)), .combine = 'c') %dopar% {
    lims <- range.chr[i, ]

    X.chr <- sub.big.matrix(X.desc,
                            firstCol = lims[1],
                            lastCol = lims[2],
                            backingpath = PATH)

    tmp <- fun.stats(X.chr, ind.train)
    S.chr <- tmp$S
    ind.col.chr <- which(tmp$pS < 10^(-thr.pvalue))
    rm(tmp)

    ind.keep <- list()
    l <- Inf
    while (l > 0) {
      ind <- ind.col.chr[which.max(S.chr[ind.col.chr])]
      ind.keep[[length(ind.keep) + 1]] <- ind

      ind.col.chr.tmp <- intersect(ind.col.chr, ind + -size:size)

      res <- R_squared_chr(pBigMat = X.chr@address,
                           rowInd = ind.train,
                           colInd = ind.col.chr.tmp,
                           colMat0 = X.chr[, ind])

      ind.del <- ind.col.chr.tmp[res > thr.corr]
      ind.col.chr <- setdiff(ind.col.chr, ind.del)
      l <- length(ind.col.chr)
      #print(l)
    }

    sort(seq2(lims)[unlist(ind.keep)])
  }
  if (!is.seq) parallel::stopCluster(cl)

  res
}

#' Title
#'
#' @param x
#' @param ind.train
#' @param size
#' @param thr.corr
#' @param ncores
#'
#' @return
#' @export
#'
#' @examples
PrunePlink <- function(x,
                  ind.train = seq(nrow(X)),
                  size = 50,
                  thr.corr = 0.5,
                  ncores = 1) {
  #check_x(x, check.y = TRUE)

  # get descriptors
  X <- x$genotypes
  X.desc <- describe(X)

  # export backingpath (force eval of promise)
  PATH <- x$backingpath

  range.chr <- LimsChr(x)

  if (is.seq <- (ncores == 1)) {
    registerDoSEQ()
  } else {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
  }
  res <- foreach(i = seq_len(nrow(range.chr)), .combine = 'c') %dopar% {
    lims <- range.chr[i, ]

    X.chr <- sub.big.matrix(X.desc,
                            firstCol = lims[1],
                            lastCol = lims[2],
                            backingpath = PATH)

    stats <- big_colstats(X.chr, ind.train)
    m.chr <- ncol(X.chr)
    keep <- rep(TRUE, m.chr)
    n <- length(ind.train)
    p <- stats$sum / (2 * n)
    maf <- pmin(p, 1 - p)
    denoX <- (n - 1) * stats$var

    keep <- R_squared_chr2(X.chr@address,
                           rowInd = ind.train,
                           keep = keep,
                           mafX = maf,
                           sumX = stats$sum,
                           denoX = denoX,
                           size = min(size, m.chr),
                           thr = thr.corr)

    seq2(lims)[which(keep)]
  }
  if (!is.seq) parallel::stopCluster(cl)

  res
}
