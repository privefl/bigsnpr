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

  stats <- bigsnpr:::bed_stats(obj.bed, ind.row, ind.col)
  af <- stats$sum / (2 * stats$nb_nona_col)
  center <- 2 * af
  scale <- sqrt(2 * af * (1 - af))

  n <- length(ind.row)
  m <- length(ind.col)

  scale2 <- scale * sqrt(stats$nb_nona_col / n)
  scale3 <- sqrt(stats$nb_nona_row / m)

  K <- FBM(n, n, init = 0)

  intervals <- bigsnpr:::CutBySize(m, block.size)

  for (j in rows_along(intervals)) {
    ind <- bigsnpr:::seq2(intervals[j, ])
    tmp <- read_bed_scaled(obj.bed, ind.row, ind.col[ind],
                           center[ind], scale2[ind])
    bigstatsr:::incrSup2(K, tcrossprod(tmp))
  }

  # Complete the lower part of the symmetric matrix
  bigstatsr:::complete2(K)
  bigstatsr:::scaleK(K, rep(0, n), rep(0, n), scale3, 0)
  structure(K, stats = stats)
}
