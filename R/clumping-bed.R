################################################################################

#' @rdname snp_clumping
#'
#' @export
#'
bed_clumping <- function(obj.bed,
                         ind.row = rows_along(obj.bed),
                         S = NULL,
                         thr.r2 = 0.2,
                         size = 100 / thr.r2,
                         exclude = NULL,
                         ncores = 1) {

  infos.chr <- obj.bed$map$chromosome
  infos.pos <- obj.bed$map$physical.pos

  check_args()

  args <- as.list(environment())

  if (!is.null(S)) assert_lengths(infos.chr, S)

  do.call(what = snp_split, args = c(args, FUN = bedClumpingChr, combine = 'c'))
}

################################################################################

bedClumpingChr <- function(obj.bed, S, ind.chr, ind.row, size, infos.pos,
                           thr.r2, exclude) {

  ind.chr <- setdiff(ind.chr, exclude)

  # cache some computations
  stats <- bed_stats(obj.bed, ind.row, ind.chr)
  center <- stats$sum / stats$nb_nona_col
  scale <- sqrt(stats$var * (stats$nb_nona_col - 1))

  # statistic to prioritize SNPs
  if (is.null(S)) {
    # use minor allele count (MAC) by default
    S.chr <- pmin(stats$sum, 2 * length(ind.row) - stats$sum)
  } else {
    S.chr <- S[ind.chr]
  }

  pos.chr <- infos.pos[ind.chr]
  assert_sorted(pos.chr)

  # main algo
  keep <- bed_clumping_chr(
    obj_bed = obj.bed,
    ind_row = ind.row,
    ind_col = ind.chr,
    center  = center,
    scale   = scale,
    ordInd  = order(S.chr, decreasing = TRUE),
    pos     = pos.chr,
    size    = size * 1000L, # in bp
    thr     = thr.r2
  )

  ind.chr[keep]
}

################################################################################
