#################################################################################

#' LD pruning and clumping
#'
#' For a `bigSNP`:
#' - `snp_pruning`: LD pruning. Similar to `--indep-pairwise size 1 thr.corr` in
#' [PLINK](http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#prune)
#' (`step` is fixed to 1).
#' - `snp_pruning`: LD clumping (and thresholding).
#' - `snp_indLRLDR`: Get SNP indices of long-range LD regions.
#'
#' @inheritParams bigsnpr-package
#'
#' @param fun.stats A function that takes a `big.matrix` __`X`__ and
#' __`ind.train`__ as parameters and returns a named list of __`S`__ and
#' __`pS`__ for every column, which are statistics and associated p-values.
#' **__`pS`__ must be computed through a decreasing function of `S`**.
#' For example, if __`S`__ follows the standard normal distribution,
#' you should use `abs(S)` instead and compute
#' `pS = 2 * pnorm(abs(S), lower.tail = FALSE)`.
#'
#' @param size
#' This parameter should be adjusted with respect to the number of SNPs.
#' - for clumping: __Radius__ of the window's size for the LD evaluations.
#' Default is `500` (I use this for a chip of 500K SNPs).
#' - for pruning: __Diameter__ of the window's size for the LD evaluations.
#' Default is `50` (as in PLINK).
#'
#'
#' @param thr.pvalue Threshold on \eqn{-log_{10}(p-value)} to assess
#' which SNPs are kept. A value of `1` is very conservative and
#' has the purpose to accelerate computations.
#' The default value of `0` means no thresholding.
#'
#' @param thr.corr Threshold on the correlation between two SNPs.
#' Default is `0.2`.
#'
#' @param exclude Vector of indices of SNPs to exclude anyway. For example,
#' can be used to exclude long-range LD regions (see Price2008).
#'
#' @references Price AL, Weale ME, Patterson N, et al.
#' Long-Range LD Can Confound Genome Scans in Admixed Populations.
#' Am J Hum Genet. 2008;83(1):132-135.
#' \url{http://dx.doi.org/10.1016/j.ajhg.2008.06.005}.
#'
#' @references Purcell, S., Neale, B., Todd-Brown, K., Thomas, L., Ferreira,
#' M., & Bender, D. et al. (2007). PLINK: A Tool Set for
#' Whole-Genome Association and Population-Based Linkage Analyses.
#' The American Journal Of Human Genetics, 81(3), 559-575.
#' \url{http://dx.doi.org/10.1086/519795}.
#'
#' @details I recommend to use clumping rather than pruning. See
#' \url{https://privefl.github.io/bigsnpr/articles/pruning-vs-clumping.html}.
#'
#' @name pruning-clumping
NULL


#' @export
#' @rdname pruning-clumping
snp_clumping <- function(x,
                         fun.stats,
                         ind.train = seq(nrow(X)),
                         thr.pvalue = 0,
                         size = 500,
                         thr.corr = 0.2,
                         exclude = NULL,
                         ncores = 1) {
  check_x(x)

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
  res <- foreach(ic = seq_len(nrow(range.chr)), .combine = 'c') %dopar% {
    lims <- range.chr[ic, ]

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

#' @export
#' @rdname pruning-clumping
snp_clumping2 <- function(x, S,
                          ind.train = seq(nrow(X)),
                          size = 1000,
                          thr.corr = 0.2,
                          exclude = NULL,
                          ncores = 1) {
  check_x(x)

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
  res <- foreach(ic = seq_len(nrow(range.chr)), .combine = 'c') %dopar% {
    lims <- range.chr[ic, ]

    X.chr <- sub.big.matrix(X.desc,
                            firstCol = lims[1],
                            lastCol = lims[2],
                            backingpath = PATH)

    # init
    ind.chr <- seq2(lims)
    ord.chr <- order(S[ind.chr], decreasing = TRUE)
    m.chr <- length(ind.chr)
    remain <- rep(TRUE, m.chr)
    remain[match(exclude, ind.chr)] <- FALSE

    # cache some computations
    stats <- bigstatsr::big_colstats(X.chr, ind.train)
    n <- length(ind.train)
    p <- stats$sum / (2 * n)
    denoX <- (n - 1) * stats$var

    # s <- setdiff(-size:size, 0)
    # keep <- rep(FALSE, m.chr)
    # for (ind in ord.chr) {
    #   if (remain[ind]) { # already excluded?
    #     ind.col <- intersect(ind + s, which(remain))
    #
    #     res <- R_squared_chr(pBigMat = X.chr@address,
    #                          rowInd = ind.train,
    #                          colInd = ind.col,
    #                          colMat0 = X.chr[, ind])
    #
    #     remain[ind.col[res > thr.corr]] <- FALSE
    #     keep[ind] <- TRUE
    #     remain[ind] <- FALSE
    #   }
    # }

    keep <- R_squared_chr3(X.chr@address,
                           rowInd = ind.train,
                           colInd = ord.chr,
                           remain = remain,
                           sumX = stats$sum,
                           denoX = denoX,
                           size = size,
                           thr = thr.corr)

    ind.chr[keep]
  }
  if (!is.seq) parallel::stopCluster(cl)

  res
}

################################################################################

#' @export
#' @rdname pruning-clumping
snp_pruning <- function(x,
                        ind.train = seq(nrow(X)),
                        size = 50,
                        thr.corr = 0.2,
                        exclude = NULL,
                        ncores = 1) {
  check_x(x)

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
  res <- foreach(ic = seq_len(nrow(range.chr)), .combine = 'c') %dopar% {
    lims <- range.chr[ic, ]

    X.chr <- sub.big.matrix(X.desc,
                            firstCol = lims[1],
                            lastCol = lims[2],
                            backingpath = PATH)

    stats <- bigstatsr::big_colstats(X.chr, ind.train)
    ind.chr <- seq2(lims)
    m.chr <- length(ind.chr)
    keep <- rep(TRUE, m.chr)
    keep[match(exclude, ind.chr)] <- FALSE
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

    ind.chr[keep]
  }
  if (!is.seq) parallel::stopCluster(cl)

  res
}

################################################################################

#' @param LD.regions A `data.frame` with columns "Chr", "Start" and "Stop".
#' Default use the table of 34 long-range LD regions that you can find
#' [there](https://goo.gl/0Ou7uI).
#'
#' @return A vector of SNP indices to be used as the __`exclude`__ parameter of
#' `snp_pruning` or `snp_clumping`.
#'
#' @import foreach
#' @export
#' @rdname pruning-clumping
snp_indLRLDR <- function(x, LD.regions = LD.wiki34) {
  chrs <- x$map$chromosome
  pos <- x$map$physical.pos

  foreach(i = 1:nrow(LD.regions), .combine = 'c') %do% {
    which((chrs == LD.regions[i, "Chr"]) &
            (pos >= LD.regions[i, "Start"]) &
            (pos <= LD.regions[i, "Stop"]))
  }
}

################################################################################
