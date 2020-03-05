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
  if (!is.null(S)) assert_lengths(infos.chr, S)

  res <- list()
  ind.chrs <- split(seq_along(infos.chr), infos.chr)
  for (ic in seq_along(ind.chrs)) {
    res[[ic]] <- bedClumpingChr(obj.bed, S, ind.chrs[[ic]], ind.row,
                                size, infos.pos, thr.r2, exclude, ncores)
  }

  sort(unlist(res))
}

################################################################################

bedClumpingChr <- function(obj.bed, S, ind.chr, ind.row, size, infos.pos,
                           thr.r2, exclude, ncores) {

  ind.chr <- setdiff(ind.chr, exclude)
  if (length(ind.chr) == 0) return(integer(0))

  # cache some computations
  stats <- bed_stats(obj.bed, ind.row, ind.chr)
  center <- stats$sum / stats$nb_nona_col
  scale <- sqrt(stats$var * (stats$nb_nona_col - 1))

  # statistic to prioritize SNPs
  if (is.null(S)) {
    # use minor allele count (MAC) by default
    S.chr <- pmin(stats$sum, 2 * stats$nb_nona_col - stats$sum)
  } else {
    S.chr <- S[ind.chr]
  }

  pos.chr <- infos.pos[ind.chr]
  assert_sorted(pos.chr)

  keep <- FBM(1, length(ind.chr), type = "integer", init = -1)
  ord <- order(S.chr, decreasing = TRUE)

  # main algo
  bed_clumping_chr(
    obj.bed,
    keep,
    ind_row = ind.row,
    ind_col = ind.chr,
    center  = center,
    scale   = scale,
    ordInd  = ord,
    rankInd = match(seq_along(ord), ord),
    pos     = pos.chr,
    size    = size * 1000L, # in bp
    thr     = thr.r2,
    ncores  = ncores
  )

  stopifnot(all(keep[] %in% 0:1))

  ind.chr[keep[] == 1]
}

################################################################################
