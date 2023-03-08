################################################################################

# Transform negative or boolean indices to positive indices
transform_ind <- function(k, lim) {

  if (missing(k))
    return(seq_len(lim))

  if (is.character(k))
    stop2("Character subsetting is not allowed.")

  res <- seq_len(lim)[k]

  if (anyNA(res))
    stop2("Error when subsetting (missing values? out of bounds?)")

  res
}


bed_accessor <- function(x, i, j, ..., drop = TRUE) {

  bigassertr::assert_nodots()

  nargs <- nargs() - !missing(drop)

  if (nargs == 2) {  # only i
    if (missing(i)) {
      nargs <- 3  # x[] is the same as x[,]
    } else {
      stop2("Vector subsetting is not allowed.")
    }
  }

  res <- read_bed(x, transform_ind(i, nrow(x)), transform_ind(j, ncol(x)))

  `if`(drop, drop(res), res)
}

################################################################################

#' Accessor methods for class `bed`.
#'
#' @param x A [bed][bed-class] object.
#' @param i A vector of indices (or nothing). You can use positive and negative
#'   indices, and also logical indices (that are recycled).
#' @param j A vector of indices (or nothing). You can use positive and negative
#'   indices, and also logical indices (that are recycled).
#' @param ... Not used. Just to make [nargs] work.
#' @param drop Whether to drop dimensions (to a vector) when a dimension is 1.
#'   Default is `TRUE`.
#'
#' @export
#'
#' @examples
#' bedfile <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
#' (obj.bed <- bed(bedfile))
#' obj.bed[1:5, 1]
#' obj.bed[1:5, 1:2]
#' typeof(obj.bed[1, 1])
#'
setMethod('[', signature(x = "bed"), bed_accessor)

################################################################################
