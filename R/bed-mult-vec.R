################################################################################

#' @description Cross-product between a "bed" object and a vector.
#'
#' Missing values are replaced by 0 (after centering), as if they
#'   had been imputed using parameter `center`.
#'
#' @param X A [bed] object.
#' @inherit bigstatsr::big_cprodVec.FBM title params return
#'
#' @export
#'
#' @examples
#' bedfile <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
#' obj.bed <- bed(bedfile)
#'
#' y.row <- rep(1, nrow(obj.bed))
#' big_cprodVec(obj.bed, y.row)
#'
big_cprodVec.bed <- function(X, y.row,
                             ind.row = rows_along(X),
                             ind.col = cols_along(X),
                             center = NULL,
                             scale = NULL) {

  assert_lengths(y.row, ind.row)

  center <- `if`(is.null(center), rep(0, length(ind.col)), as_vec(center))
  assert_lengths(center, ind.col)

  scale <- `if`(is.null(scale), rep(1, length(ind.col)), as_vec(scale))
  assert_lengths(scale, ind.col)

  bed_cpMatVec4(X, ind.row, ind.col, center, scale, y.row)
}

################################################################################

#' @description Product between a "bed" object and a vector.
#'
#' Missing values are replaced by 0 (after centering), as if they
#'   had been imputed using parameter `center`.
#'
#' @param X A [bed] object.
#' @inherit bigstatsr::big_prodVec.FBM title params return
#'
#' @export
#'
#' @examples
#' bedfile <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
#' obj.bed <- bed(bedfile)
#'
#' y.col <- rep(1, ncol(obj.bed))
#' big_prodVec(obj.bed, y.col)
#'
big_prodVec.bed <- function(X, y.col,
                            ind.row = rows_along(X),
                            ind.col = cols_along(X),
                            center = NULL,
                            scale = NULL) {

  assert_lengths(y.col, ind.col)

  center <- `if`(is.null(center), rep(0, length(ind.col)), as_vec(center))
  assert_lengths(center, ind.col)

  scale <- `if`(is.null(scale), rep(1, length(ind.col)), as_vec(scale))
  assert_lengths(scale, ind.col)

  bed_pMatVec4(X, ind.row, ind.col, center, scale, y.col)
}

################################################################################
