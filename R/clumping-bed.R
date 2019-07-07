################################################################################

#' @rdname snp_clumping
#'
#' @export
#'
bed_clumping <- function(bedfile,
                         ind.row = rows_along(obj.bed),
                         S = NULL,
                         thr.r2 = 0.2,
                         size = 100 / thr.r2,
                         exclude = NULL,
                         ncores = 1) {

  obj.bed <- bed(bedfile)
  rm(bedfile)

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
    af <- center / 2
    S.chr <- pmin(af, 1 - af)
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

#' @import foreach
#' @export
#' @rdname snp_clumping
bed_indLRLDR <- function(bedfile, LD.regions = LD.wiki34) {

  snp_infos <- bigreadr::fread2(sub_bed(bedfile, ".bim"))
  infos.chr <- snp_infos[[1]]
  infos.pos <- snp_infos[[4]]

  check_args()

  foreach(ic = rows_along(LD.regions), .combine = 'c') %do% {
    which((infos.chr == LD.regions[ic, "Chr"]) &
            (infos.pos >= LD.regions[ic, "Start"]) &
            (infos.pos <= LD.regions[ic, "Stop"]))
  }
}

################################################################################
