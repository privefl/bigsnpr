################################################################################

#' @inherit bigstatsr::big_tcrossprodSelf title description params return details
#' @inheritParams bigsnpr-package
#'
#' @export
#'
#' @examples
#' bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
#' obj.bed <- bed(bedfile)
#'
#' K <- bed_tcrossprodSelf(obj.bed)
#' K[1:4, 1:6] / ncol(obj.bed)
#'
bed_tcrossprodSelf <- function(
  obj.bed,
  fun.scaling = bed_scaleBinom,
  ind.row = rows_along(obj.bed),
  ind.col = cols_along(obj.bed),
  block.size = block_size(length(ind.row))
) {

  n <- length(ind.row)
  K <- FBM(n, n, init = 0)

  m <- length(ind.col)
  center <- numeric(m)
  scale  <- numeric(m)

  intervals <- CutBySize(m, block.size)

  for (j in rows_along(intervals)) {
    ind <- seq2(intervals[j, ])
    ind.col.ind <- ind.col[ind]
    ms <- fun.scaling(obj.bed, ind.row = ind.row, ind.col = ind.col.ind)
    center[ind] <- ms$center
    scale[ind]  <- ms$scale
    tmp <- read_bed_scaled(obj.bed, ind.row, ind.col.ind,
                           ms$center, ms$scale)
    big_increment(K, tcrossprod(tmp))
  }

  structure(K, center = center, scale = scale)
}

################################################################################
