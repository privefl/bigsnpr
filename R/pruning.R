#################################################################################

#' LD pruning and clumping
#'
#' For a `bigSNP`:
#' - `snp_pruning`: LD pruning. Similar to "`--indep-pairwise size 1 thr.corr`"
#' in [PLINK](http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#prune)
#' (`step` is fixed to 1).
#' - `snp_clumping`: LD clumping.
#' - `snp_indLRLDR`: Get SNP indices of long-range LD regions.
#'
#' @inheritParams bigsnpr-package
#'
#' @param S A vector of column statistics which express the importance
#' of each SNP (the more important is the SNP, the greater should be
#' the corresponding statistic).
#' For example, if `S` follows the standard normal distribution,
#' you should probably use `abs(S)` instead.
#'
#' @param size
#' This parameter should be adjusted with respect to the number of SNPs.
#' - for clumping: __Radius__ of the window's size for the LD evaluations.
#' Default is `500` (I use this for a chip of 500K SNPs).
#' - for pruning: __Diameter__ of the window's size for the LD evaluations.
#' Default is `50` (as in PLINK).
#'
#' @param thr.corr Threshold on the correlation between two SNPs.
#' Default is `0.5`.
#'
#' @param exclude Vector of indices of SNPs to exclude anyway. For example,
#' can be used to exclude long-range LD regions (see Price2008). Another use
#' can be for thresholding with respect to p-values associated with `S`.
#'
#' @references Price AL, Weale ME, Patterson N, et al.
#' Long-Range LD Can Confound Genome Scans in Admixed Populations.
#' Am J Hum Genet. 2008;83(1):132-135.
#' [DOI](http://dx.doi.org/10.1016/j.ajhg.2008.06.005)
#'
#' @return
#' - `snp_pruning` & `snp_pruning`: SNP indices which are __kept__.
#' - `snp_indLRLDR`: SNP indices to be used as (part of) the
#' __`exclude`__ parameter of `snp_pruning` or `snp_clumping`.
#'
#' @details I recommend to use clumping rather than pruning. See
#' [this article](https://privefl.github.io/bigsnpr/articles/pruning-vs-clumping.html).
#'
#' @name pruning-clumping
NULL

#################################################################################

#' @export
#' @rdname pruning-clumping
snp_clumping <- function(x, S,
                         ind.train = seq(nrow(X)),
                         size = 500,
                         thr.corr = 0.5,
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
    remain <- rep(TRUE, length(ind.chr))
    remain[match(exclude, ind.chr)] <- FALSE

    # cache some computations
    stats <- bigstatsr::big_colstats(X.chr, ind.train)
    n <- length(ind.train)
    p <- stats$sum / (2 * n)
    denoX <- (n - 1) * stats$var


    keep <- clumping(X.chr@address,
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
                        thr.corr = 0.5,
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

    keep <- pruning(X.chr@address,
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
