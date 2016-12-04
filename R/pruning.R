#################################################################################

#' LD clumping
#'
#' LD clumping (and thresholding) for a `bigSNP`.
#'
#' @inheritParams bigsnpr-package
#'
#' @param fun.stats A function that takes a `big.matrix` __`X`__ and
#' __`ind.train`__ as parameters and returns a named list of __`S`__ and
#' __`pS`__ for every column, which are statistics and associated p-values.
#' **__`pS`__ needs to be computed through a decreasing function of `S`**.
#' For example, if __`S`__ follows the standard normal distribution,
#' you should use `abs(S)` instead and compute
#' `pS = 2 * pnorm(abs(S), lower.tail = FALSE)`.
#'
#' @param size __Radius__ of the window's size for the LD evaluations.
#' Default is `500` for clumping. This should be
#' adjusted for different number of SNPs (this corresponds to the defaults
#' I use for a chip of 500K SNPs).
#'
#' @param thr.pvalue Threshold on \eqn{-log_{10}(p-value)} to assess
#' which SNPs are kept. A value of `1` is very conservative and
#' has the purpose to accelerate computations.
#' A value of `0` means no thresholding and is the default.
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
#' @details I recommend to use clumping rather than pruning. See
#' \url{https://privefl.github.io/bigsnpr/articles/pruning-vs-clumping.html}.
#'
#' @export
snp_clumping <- function(x,
                  ind.train = seq(nrow(X)),
                  fun.stats,
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

################################################################################

#' LD pruning
#'
#' LD pruning (and thresholding) for a `bigSNP`.
#' Similar to `--indep-pairwise size 1 thr.corr` (`step` is fixed to 1).
#'
#' @inherit snp_clumping params references
#' @param size __Diameter__ of the window's size for the LD evaluations.
#' Default is `50` for pruning (default in PLINK) and should be adjusted
#' for different number of SNPs.
#'
#' @references Purcell, S., Neale, B., Todd-Brown, K., Thomas, L., Ferreira,
#' M., & Bender, D. et al. (2007). PLINK: A Tool Set for
#' Whole-Genome Association and Population-Based Linkage Analyses.
#' The American Journal Of Human Genetics, 81(3), 559-575.
#' \url{http://dx.doi.org/10.1086/519795}.
#'
#' @export
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

    stats <- big_colstats(X.chr, ind.train)
    m.chr <- ncol(X.chr)
    keep <- rep(TRUE, m.chr)
    keep[match(exclude, seq2(lims))] <- FALSE
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

################################################################################

#' Get SNPs of long-range LD regions
#'
#' @inheritParams bigsnpr-package
#'
#' @param LD.regions A `data.frame` with columns "Chr", "Start" and "Stop".
#' Default use the table of 34 long-range LD regions that you can find there:
#' \url{https://goo.gl/0Ou7uI}.
#'
#' @return A vector of SNP indices to be used as the __`exclude`__ parameter of
#' `snp_pruning` or `snp_clumping`.
#' @export
#' @import foreach
#'
#' @examples
excludeLDreg <- function(x, LD.regions = LD.wiki34) {
  chrs <- x$map$chromosome
  pos <- x$map$physical.pos

  foreach(i = 1:nrow(LD.regions), .combine = 'c') %do% {
    chr = LD.regions[i, "Chr"]
    start = LD.regions[i, "Start"]
    end = LD.regions[i, "Stop"]
    which((chrs == chr) & (pos >= start) & (pos <= end))
  }
}

################################################################################
