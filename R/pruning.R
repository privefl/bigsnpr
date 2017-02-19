#################################################################################

#' LD pruning and clumping
#'
#' For a `bigSNP`:
#' - `snp_pruning`: LD pruning. Similar to "`--indep-pairwise size 1 thr.r2`" in
#'   [PLINK 1.07](http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#prune)
#'   (`step` is fixed to 1).
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
#' Default is `1000` (I use this for a chip of 500K SNPs).
#' - for pruning: __Diameter__ of the window's size for the LD evaluations.
#' Default is `50` (as in PLINK 1.07).
#'
#' @param thr.r2 Threshold over the squared correlation between two SNPs.
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
#' @example examples/example-pruning.R
#'
#' @name pruning-clumping
NULL

#################################################################################

#' @export
#' @rdname pruning-clumping
snp_clumping <- function(x, S,
                         ind.train = seq(nrow(X)),
                         size = 1000,
                         thr.r2 = 0.5,
                         exclude = NULL,
                         ncores = 1) {
  check_x(x)

  # get descriptors
  X <- x$genotypes
  X.desc <- describe(X)

  # get ranges of chromosomes
  range.chr <- LimsChr(x)

  if (is.seq <- (ncores == 1)) {
    registerDoSEQ()
  } else {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
  }
  res <- foreach(ic = seq_len(nrow(range.chr)), .combine = 'c') %dopar% {
    lims <- range.chr[ic, ]

    X <- attach.big.matrix(X.desc)

    # init
    ind.chr <- seq2(lims)
    ord.chr <- order(S[ind.chr], decreasing = TRUE)
    remain <- rep(TRUE, length(ind.chr))
    remain[match(exclude, ind.chr)] <- FALSE

    # cache some computations
    stats <- bigstatsr::big_colstats(X, ind.train, ind.chr)
    n <- length(ind.train)
    denoX <- (n - 1) * stats$var

    # main algo
    keep <- clumping(X@address,
                     rowInd = ind.train,
                     colInd = ind.chr,
                     ordInd = ord.chr,
                     remain = remain,
                     sumX = stats$sum,
                     denoX = denoX,
                     size = size,
                     thr = thr.r2)

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
                        is.size.in.kb = FALSE,
                        thr.r2 = 0.5,
                        exclude = NULL,
                        ncores = 1) {
  check_x(x)

  # get descriptors
  X <- x$genotypes
  X.desc <- describe(X)

  # get ranges of chromosomes
  range.chr <- LimsChr(x)
  pos <- x$map$physical.pos

  if (is.seq <- (ncores == 1)) {
    registerDoSEQ()
  } else {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }
  res <- foreach(ic = seq_len(nrow(range.chr)), .combine = 'c') %dopar% {
    lims <- range.chr[ic, ]

    X <- attach.big.matrix(X.desc)

    # cache some computations
    ind.chr <- seq2(lims)
    stats <- bigstatsr::big_colstats(X, ind.train, ind.chr)
    m.chr <- length(ind.chr)
    keep <- rep(TRUE, m.chr)
    keep[match(exclude, ind.chr)] <- FALSE
    n <- length(ind.train)
    p <- stats$sum / (2 * n)
    maf <- pmin(p, 1 - p)
    denoX <- (n - 1) * stats$var
    nulls <- which(denoX == 0)
    if (l <- length(nulls)) {
      message(sprintf("Excluding %d monoallelic markers...", l))
      keep[nulls] <- FALSE
    }


    # main algo
    if (is.size.in.kb) {
      keep <- pruning2(X@address,
                       rowInd = ind.train,
                       colInd = ind.chr,
                       keep = keep,
                       pos = c(pos[ind.chr], .Machine$integer.max),
                       mafX = maf,
                       sumX = stats$sum,
                       denoX = denoX,
                       size = size * 1000, # in kb
                       thr = thr.r2)
    } else {
      keep <- pruning(X@address,
                      rowInd = ind.train,
                      colInd = ind.chr,
                      keep = keep,
                      mafX = maf,
                      sumX = stats$sum,
                      denoX = denoX,
                      size = min(size, m.chr),
                      thr = thr.r2)
    }


    ind.chr[keep]
  }
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
