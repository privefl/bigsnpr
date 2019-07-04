################################################################################

#' LD clumping
#'
#' For a `bigSNP`:
#' - `snp_pruning()`: LD pruning. Similar to "`--indep-pairwise (size+1) 1 thr.r2`"
#'   in [PLINK](https://www.cog-genomics.org/plink/1.9/ld)
#'   (`step` is fixed to 1). **This function is deprecated (see
#' [this article](https://privefl.github.io/bigsnpr/articles/pruning-vs-clumping.html)).**
#' - `snp_clumping()` (and `bed_clumping()`): LD clumping. If you do not provide any statistic to rank
#'   SNPs, it would use minor allele frequencies (MAFs), making clumping similar
#'   to pruning.
#' - `snp_indLRLDR()`: Get SNP indices of long-range LD regions for the
#'   human genome.
#'
#' @inheritParams bigsnpr-package
#' @inheritParams snp_readBed
#'
#' @param S A vector of column statistics which express the importance
#' of each SNP (the more important is the SNP, the greater should be
#' the corresponding statistic).\cr
#' For example, if `S` follows the standard normal distribution, and "important"
#' means significantly different from 0, you must use `abs(S)` instead.\cr
#' If not specified, the MAF is computed and used.
#'
#' @param size For one SNP, window size around this SNP to compute correlations.
#' Default is `100 / thr.r2` for clumping (0.2 -> 500; 0.1 -> 1000; 0.5 -> 200).
#' If not providing `infos.pos` (`NULL`, the default), this is a window in
#' number of SNPs, otherwise it is a window in kb (genetic distance).
#' I recommend that you provide the positions if available.
#'
#' @param exclude Vector of SNP indices to exclude anyway. For example,
#' can be used to exclude long-range LD regions (see Price2008). Another use
#' can be for thresholding with respect to p-values associated with `S`.
#'
#' @param thr.r2 Threshold over the squared correlation between two SNPs.
#'   Default is `0.2`.
#'
#' @references Price AL, Weale ME, Patterson N, et al.
#' Long-Range LD Can Confound Genome Scans in Admixed Populations.
#' Am J Hum Genet. 2008;83(1):132-135.
#' \url{http://dx.doi.org/10.1016/j.ajhg.2008.06.005}
#'
#' @return
#' - `snp_clumping()` (and `bed_clumping()`): SNP indices that are __kept__.
#' - `snp_indLRLDR()`: SNP indices to be used as (part of) the '__`exclude`__'
#'   parameter of `snp_clumping()`.
#'
#' @example examples/example-pruning.R
#'
#' @export
#'
snp_clumping <- function(G, infos.chr,
                         ind.row = rows_along(G),
                         S = NULL,
                         thr.r2 = 0.2,
                         size = 100 / thr.r2,
                         infos.pos = NULL,
                         is.size.in.bp = NULL,
                         exclude = NULL,
                         ncores = 1) {

  check_args()

  if (!missing(is.size.in.bp))
    warning2("Parameter 'is.size.in.bp' is deprecated.")

  args <- as.list(environment())

  if (!is.null(S)) assert_lengths(infos.chr, S)

  do.call(what = snp_split, args = c(args, FUN = clumpingChr, combine = 'c'))
}

################################################################################

clumpingChr <- function(G, S, ind.chr, ind.row, size, is.size.in.bp, infos.pos,
                        thr.r2, exclude) {

  ind.chr <- setdiff(ind.chr, exclude)

  # cache some computations
  stats <- big_colstats(G, ind.row = ind.row, ind.col = ind.chr)
  n <- length(ind.row)

  # statistic to prioritize SNPs
  if (is.null(S)) {
    af <- stats$sum / (2 * n)
    S.chr <- pmin(af, 1 - af)
  } else {
    S.chr <- S[ind.chr]
  }
  pos.chr <- `if`(is.null(infos.pos), 1000L * seq_along(ind.chr), infos.pos[ind.chr])
  assert_sorted(pos.chr)

  # main algo
  keep <- clumping_chr(
    G,
    rowInd = ind.row,
    colInd = ind.chr,
    ordInd = order(S.chr, decreasing = TRUE),
    pos    = pos.chr,
    sumX   = stats$sum,
    denoX  = (n - 1) * stats$var,
    size   = size * 1000L, # in bp
    thr    = thr.r2
  )

  ind.chr[keep]
}

################################################################################

#' @export
#' @rdname snp_clumping
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

  warning2("Pruning is deprecated; using clumping (on MAF) instead..\n%s",
           "See why there: https://goo.gl/Td5YYv.")

  args <- c(as.list(environment()), list(S = NULL))
  args[["nploidy"]] <- NULL

  do.call(what = snp_split, args = c(args, FUN = clumpingChr, combine = 'c'))
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
#' @rdname snp_clumping
snp_indLRLDR <- function(infos.chr, infos.pos, LD.regions = LD.wiki34) {

  check_args()

  foreach(ic = rows_along(LD.regions), .combine = 'c') %do% {
    which((infos.chr == LD.regions[ic, "Chr"]) &
            (infos.pos >= LD.regions[ic, "Start"]) &
            (infos.pos <= LD.regions[ic, "Stop"]))
  }
}

################################################################################
