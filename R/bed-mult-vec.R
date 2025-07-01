################################################################################

#' @description Cross-product between a "bed" object and a vector.
#'
#' Missing values are replaced by 0 (after centering), as if they
#'   had been imputed using parameter `center`.
#'
#' @inheritParams bigsnpr-package
#' @inherit bigstatsr::big_cprodVec title params return
#'
#' @export
#'
#' @examples
#' bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
#' obj.bed <- bed(bedfile)
#'
#' y.row <- rep(1, nrow(obj.bed))
#' str(bed_cprodVec(obj.bed, y.row))
#'
bed_cprodVec <- function(obj.bed, y.row,
                         ind.row = rows_along(obj.bed),
                         ind.col = cols_along(obj.bed),
                         center = rep(0, length(ind.col)),
                         scale  = rep(1, length(ind.col)),
                         ncores = 1) {

  check_args()
  assert_lengths(y.row, ind.row)

  center <- as_vec(center)
  assert_lengths(center, ind.col)

  scale <- as_vec(scale)
  assert_lengths(scale, ind.col)

  bed_cpMatVec4(obj.bed, ind.row, ind.col, center, scale, y.row, ncores)
}

################################################################################

#' @description Product between a "bed" object and a vector.
#'
#' Missing values are replaced by 0 (after centering), as if they
#'   had been imputed using parameter `center`.
#'
#' @inheritParams bigsnpr-package
#' @inherit bigstatsr::big_prodVec title params return
#'
#' @export
#'
#' @examples
#' bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
#' obj.bed <- bed(bedfile)
#'
#' y.col <- rep(1, ncol(obj.bed))
#' str(bed_prodVec(obj.bed, y.col))
#'
bed_prodVec <- function(obj.bed, y.col,
                        ind.row = rows_along(obj.bed),
                        ind.col = cols_along(obj.bed),
                        center = rep(0, length(ind.col)),
                        scale = rep(1, length(ind.col)),
                        ncores = 1) {

  check_args()
  assert_lengths(y.col, ind.col)

  center <- as_vec(center)
  assert_lengths(center, ind.col)

  scale <- as_vec(scale)
  assert_lengths(scale, ind.col)

  bed_pMatVec4(obj.bed, ind.row, ind.col, center, scale, y.col, ncores)
}

################################################################################
