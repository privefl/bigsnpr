################################################################################

#' LD clumping
#'
#' For a `bigSNP`:
#' - `snp_pruning()`: LD pruning. Similar to "`--indep-pairwise (size+1) 1 thr.r2`"
#'   in [PLINK](https://www.cog-genomics.org/plink/1.9/ld).
#'   **This function is deprecated (see [this article](http://bit.ly/2uKo3MN)).**
#' - `snp_clumping()` (and `bed_clumping()`): LD clumping. If you do not provide
#'   any statistic to rank SNPs, it would use minor allele frequencies (MAFs),
#'   making clumping similar to pruning.
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
#' **If not specified, MAFs are computed and used.**
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
#' @export
#'
#' @examples
#' test <- snp_attachExtdata()
#' G <- test$genotypes
#'
#' # clumping (prioritizing higher MAF)
#' ind.keep <- snp_clumping(G, infos.chr = test$map$chromosome,
#'                          infos.pos = test$map$physical.pos,
#'                          thr.r2 = 0.1)
#'
#' # keep most of them -> not much LD in this simulated dataset
#' length(ind.keep) / ncol(G)
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

  if (!is.null(S)) assert_lengths(infos.chr, S)

  ind.noexcl <- setdiff(seq_along(infos.chr), exclude)

  sort(unlist(
    lapply(split(ind.noexcl, infos.chr), function(ind.chr) {
      clumpingChr(G, S, ind.chr, ind.row, size, infos.pos, thr.r2, ncores)
    }),
    use.names = FALSE
  ))
}

################################################################################

clumpingChr <- function(G, S, ind.chr, ind.row, size, infos.pos, thr.r2, ncores) {

  # cache some computations
  stats <- big_colstats(G, ind.row = ind.row, ind.col = ind.chr, ncores = ncores)

  # statistic to prioritize SNPs
  n <- length(ind.row)
  if (is.null(S)) {
    af <- stats$sum / (2 * n)
    S.chr <- pmin(af, 1 - af)
  } else {
    S.chr <- S[ind.chr]
  }
  ord <- order(S.chr, decreasing = TRUE)

  if (is.null(infos.pos)) {
    pos.chr <- seq_along(ind.chr)
  } else {
    size <- size * 1000  # kbp to bp
    pos.chr <- infos.pos[ind.chr]
    assert_sorted(pos.chr)
  }

  keep <- FBM(1, length(ind.chr), type = "integer", init = -1)

  # main algo
  clumping_chr(
    G,
    keep,
    rowInd = ind.row,
    colInd = ind.chr,
    ordInd = ord,
    rankInd = match(seq_along(ord), ord),
    pos    = pos.chr,
    sumX   = stats$sum,
    denoX  = (n - 1) * stats$var,
    size   = size,
    thr    = thr.r2,
    ncores = ncores
  )

  stopifnot(all(keep[] %in% 0:1))

  ind.chr[keep[] == 1]
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
                        nploidy = 2,
                        ncores = 1) {

  stop2("Pruning is deprecated; please use clumping (on MAF) instead..\n%s",
        "See why at http://bit.ly/2uKo3MN.")
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
