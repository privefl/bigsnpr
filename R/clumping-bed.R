################################################################################

#' @rdname snp_clumping
#'
#' @export
#'
bed_clumping <- function(bedfile,
                         ind.row = seq_len(n_total),
                         S = NULL,
                         thr.r2 = 0.2,
                         size = 100 / thr.r2,
                         exclude = NULL,
                         ncores = 1) {

  # Check extension of file
  assert_ext(bedfile, "bed")
  # Get other files
  bimfile <- sub("\\.bed$", ".bim", bedfile)
  famfile <- sub("\\.bed$", ".fam", bedfile)
  # Check if all three files exist
  sapply(c(bedfile, bimfile, famfile), assert_exist)

  n_total <- bigreadr::nlines(famfile)
  snp_info <- bigreadr::fread2(bimfile, select = c(1, 4))
  m_total <- nrow(snp_info)
  infos.chr <- snp_info[[1]]
  infos.pos <- snp_info[[2]]
  rm(bimfile, famfile, snp_info)

  args <- as.list(environment())

  if (!is.null(S)) assert_lengths(infos.chr, S)

  do.call(what = snp_split, args = c(args, FUN = bedClumpingChr, combine = 'c'))
}

################################################################################

bedClumpingChr <- function(bedfile, n_total, m_total,
                           S, ind.chr, ind.row, size, infos.pos,
                           thr.r2, exclude) {

  ind.chr <- setdiff(ind.chr, exclude)

  # cache some computations
  lookup_byte <- getCode()
  stats <- bedcolvars(bedfile, n_total, m_total, ind.row, ind.chr, lookup_byte)
  center <- stats$sum / stats$nona
  scale <- sqrt(stats$var * (stats$nona - 1))
  lookup_scale <- rbind(t(sapply(0:2, function(g) (g - center) / scale)), 0)

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
    bedfile, n_total, m_total,
    ind.row, ind.chr,
    lookup_byte, lookup_scale,
    ordInd = order(S.chr, decreasing = TRUE),
    pos    = pos.chr,
    size   = size * 1000L, # in bp
    thr    = thr.r2
  )

  ind.chr[keep]
}

################################################################################
