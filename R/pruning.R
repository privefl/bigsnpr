#################################################################################

#' LD pruning and clumping
#'
#' For a `bigSNP`:
#' - `snp_pruning`: LD pruning. Similar to "`--indep-pairwise (size+1) 1 thr.r2`"
#'   in [PLINK](https://www.cog-genomics.org/plink/1.9/ld).
#'   (`step` is fixed to 1).
#' - `snp_clumping`: LD clumping.
#' - `snp_indLRLDR`: Get SNP indices of long-range LD regions for the
#'   human genome.
#'
#' @inheritParams bigsnpr-package
#'
#' @param S A vector of column statistics which express the importance
#' of each SNP (the more important is the SNP, the greater should be
#' the corresponding statistic).\cr
#' For example, if `S` follows the standard normal distribution, and "important"
#' means significantly different from 0, you must use `abs(S)` instead.\cr
#' If not specified, the MAF is computed and used.
#'
#' @param size For one SNP, number of SNPs at its left and its right to
#' be tested for being correlated with this particular SNP.
#' This parameter should be adjusted with respect to the number of SNPs.
#' The default are
#' - `49` for pruning (as in PLINK),
#' - `500` for clumping (I use this for a chip of 500K SNPs).
#'
#' @param thr.r2 Threshold over the squared correlation between two SNPs.
#' Default is `0.2`.
#'
#' @param exclude Vector of SNP indices to exclude anyway. For example,
#' can be used to exclude long-range LD regions (see Price2008). Another use
#' can be for thresholding with respect to p-values associated with `S`.
#'
#' @param is.size.in.bp Is the specified size is corresponding to
#' kilo base-pairs. Default is `FALSE` (it is only SNP counts).
#'
#' @references Price AL, Weale ME, Patterson N, et al.
#' Long-Range LD Can Confound Genome Scans in Admixed Populations.
#' Am J Hum Genet. 2008;83(1):132-135.
#' \url{http://dx.doi.org/10.1016/j.ajhg.2008.06.005}
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

clumpingChr <- function(G, S, ind.chr, ind.row, size, is.size.in.bp, infos.pos,
                        thr.r2, exclude) {

  # init
  ord.chr <- order(S[ind.chr], decreasing = TRUE)
  remain <- rep(TRUE, length(ind.chr))
  remain[match(exclude, ind.chr)] <- FALSE

  # cache some computations
  stats <- big_colstats(G, ind.row = ind.row, ind.col = ind.chr)
  n <- length(ind.row)
  denoX <- (n - 1) * stats$var

  # main algo
  if (is.size.in.bp) {
    stopifnot(length(infos.pos) != 0) # TODO: better test
    keep <- clumping2(G,
                      rowInd = ind.row,
                      colInd = ind.chr,
                      ordInd = ord.chr,
                      remain = remain,
                      pos = c(infos.pos[ind.chr], .Machine$integer.max),
                      sumX = stats$sum,
                      denoX = denoX,
                      size = size * 1000, # in bp
                      thr = thr.r2)
  } else {
    keep <- clumping(G,
                     rowInd = ind.row,
                     colInd = ind.chr,
                     ordInd = ord.chr,
                     remain = remain,
                     sumX = stats$sum,
                     denoX = denoX,
                     size = size,
                     thr = thr.r2)
  }

  ind.chr[keep]
}

#' @export
#' @rdname pruning-clumping
snp_clumping <- function(G, infos.chr,
                         ind.row = rows_along(G),
                         S = NULL,
                         size = 500,
                         is.size.in.bp = FALSE,
                         infos.pos = NULL,
                         thr.r2 = 0.2,
                         exclude = NULL,
                         ncores = 1) {

  check_args()

  if (is.null(S)) S <- snp_MAF(G, ind.row = ind.row)
  args <- as.list(environment())

  assert_lengths(infos.chr, S)

  do.call(what = snp_split, args = c(args, FUN = clumpingChr, combine = 'c'))
}

################################################################################

pruningChr <- function(G, ind.chr, ind.row, nploidy,
                       size, is.size.in.bp, infos.pos, thr.r2, exclude) {

  # cache some computations
  stats <- big_colstats(G, ind.row, ind.chr)
  m.chr <- length(ind.chr)
  keep <- rep(TRUE, m.chr)
  keep[match(exclude, ind.chr)] <- FALSE
  n <- length(ind.row)
  p <- stats$sum / (nploidy * n)
  maf <- pmin(p, 1 - p)
  denoX <- (n - 1) * stats$var
  # nulls <- which(denoX == 0)
  # if (l <- length(nulls)) {
  #   message2("Excluding %d monoallelic markers...", l)
  #   keep[nulls] <- FALSE
  # }

  # main algo
  if (is.size.in.bp) {
    stopifnot(length(infos.pos) != 0) # TODO: better test
    keep <- pruning2(G,
                     rowInd = ind.row,
                     colInd = ind.chr,
                     keep = keep,
                     pos = c(infos.pos[ind.chr], .Machine$integer.max),
                     mafX = maf,
                     sumX = stats$sum,
                     denoX = denoX,
                     size = size * 1000, # in bp
                     thr = thr.r2)
  } else {
    keep <- pruning(G,
                    rowInd = ind.row,
                    colInd = ind.chr,
                    keep = keep,
                    mafX = maf,
                    sumX = stats$sum,
                    denoX = denoX,
                    size = size,
                    thr = thr.r2)
  }

  ind.chr[keep]
}

#' @export
#' @rdname pruning-clumping
snp_pruning <- function(G, infos.chr,
                        ind.row = rows_along(G),
                        size = 49,
                        is.size.in.bp = FALSE,
                        infos.pos = NULL,
                        thr.r2 = 0.2,
                        exclude = NULL,
                        nploidy = getOption("bigsnpr.nploidy"),
                        ncores = 1) {

  check_args()

  args <- as.list(environment())

  do.call(what = snp_split, args = c(args, FUN = pruningChr, combine = 'c'))
}

################################################################################

#' Long-range LD regions
#'
#' 34 long-range Linkage Disequilibrium (LD) regions for the human genome
#' based on some [wiki table](https://goo.gl/0Ou7uI).
#'
#' @format A data frame with 34 rows (regions) and 4 variables:
#' - `Chr`: region's chromosome
#' - `Start`: starting position of the region (in bp)
#' - `Stop`: stopping position of the region (in bp)
#' - `ID`: some ID of the region.
"LD.wiki34"

#' @param LD.regions A `data.frame` with columns "Chr", "Start" and "Stop".
#' Default use the table of 34 long-range LD regions that you can find
#' [there](https://goo.gl/0Ou7uI).
#'
#' @import foreach
#' @export
#' @rdname pruning-clumping
snp_indLRLDR <- function(infos.chr, infos.pos, LD.regions = LD.wiki34) {

  check_args()

  foreach(ic = 1:nrow(LD.regions), .combine = 'c') %do% {
    which((infos.chr == LD.regions[ic, "Chr"]) &
            (infos.pos >= LD.regions[ic, "Start"]) &
            (infos.pos <= LD.regions[ic, "Stop"]))
  }
}

################################################################################
