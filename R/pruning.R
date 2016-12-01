#' LD clumping
#'
#' LD clumping (and thresholding) for a `bigSNP`.
#'
#' @name pruning
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
#' @param size Radius of the window's size for the LD evaluations.
#' Default is `1000`, which seems pretty conservative for a standard
#' chip of less than a million SNPs.
#'
#' @param thr.pvalue Threshold on \eqn{-log_{10}(p-value)} to assess
#' which SNPs are kept. A default of `1` is very conservative and
#' has the purpose to accelerate computations.
#'
#' @param thr.corr Threshold on the correlation between two SNPs.
#' SNPs which are too correlated with another SNP which is more correlated
#' with the disease are pruned.
#'
#' @param exclude Vector of indices of SNPs to exclude anyway. For example,
#' can be used to exclude long-range LD regions (see Price2008)
#'
#' @references Price AL, Weale ME, Patterson N, et al.
#' Long-Range LD Can Confound Genome Scans in Admixed Populations.
#' Am J Hum Genet. 2008;83(1):132-135.
#' \link{http://dx.doi.org/10.1016/j.ajhg.2008.06.005}.
#'
#' @example examples/example.pruning.R
NULL

#' @export
#' @rdname snp_pruning
snp_clump <- function(x,
                  ind.train = seq(nrow(X)),
                  fun.stats,
                  thr.pvalue = 1,
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
#' @rdname snp_pruning
snp_prune <- function(x,
                  ind.train = seq(nrow(X)),
                  size = 50,
                  thr.corr = 0.5,
                  exclude = NULL,
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

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
excludeLDreg <- function(x) {
  url <- "http://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_%28LD%29"
  tables <- XML::readHTMLTable(url, header = TRUE, which = 1,
                          colClasses = c(rep("integer", 3), "character"),
                          stringsAsFactors = FALSE)
  names(tables) <- c("Chr", "Start", "Stop", "ID")

  chrs <- x$map$chromosome
  pos <- x$map$physical.pos

  require(foreach)
  LD <- foreach(i = 1:nrow(tables), .combine = 'c') %do% {
    chr = tables[i, "Chr"]
    start = tables[i, "Start"]
    end = tables[i, "Stop"]
    which((chrs == chr) & (pos >= start) & (pos <= end))
  }
}
