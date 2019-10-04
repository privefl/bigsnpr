################################################################################

#' @inherit bigstatsr::big_tcrossprodSelf title description params return details
#' @inheritParams bigsnpr-package
#'
#' @export
#'
#' @examples
#' bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
#' obj.bed <- bed(bedfile)
#' K <- bed_tcrossprodSelf(obj.bed)
#' K[1:4, 1:6] / ncol(obj.bed)
bed_tcrossprodSelf <- function(
  obj.bed,
  ind.row = rows_along(obj.bed),
  ind.col = cols_along(obj.bed),
  block.size = block_size(length(ind.row))
) {

  stats <- bed_stats(obj.bed, ind.row, ind.col)
  af <- stats$sum / (2 * stats$nb_nona_col)
  center <- 2 * af
  scale <- sqrt(2 * af * (1 - af))

  n <- length(ind.row)
  m <- length(ind.col)

  K <- FBM(n, n, init = 0)

  intervals <- CutBySize(m, block.size)

  for (j in rows_along(intervals)) {
    ind <- seq2(intervals[j, ])
    tmp <- read_bed_scaled(obj.bed, ind.row, ind.col[ind],
                           center[ind], scale[ind])
    big_increment(K, tcrossprod(tmp))
  }

  structure(K, stats = stats)
}

################################################################################
