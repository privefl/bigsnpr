#' @import bigmemory bigstatsr
#' @useDynLib bigsnpr
#' @importFrom Rcpp sourceCpp
#'
#' @param X A [big.matrix][bigmemory::big.matrix-class].
#' You shouldn't have missing values in your data.
#'
#' @param x A [bigSNP][bigSNP-class].
#'
#' @param has.pheno Does __`x`__ has well-specified phenotypes?
#' i.e. `x$fam$affection` is a vector with exactly two values
#' (and the value coding for cases is larger than the value coding
#' for controls). So, in particular, no missing value is allowed.
#' Default is `TRUE`.
#'
#' @param ind.train An optional vector of the row indices that are used,
#' for the training part. If not specified, all data are used.
#'
#' @param block.size Maximum number of loci read at once (for all individuals).
#'
#' @param ncores Number of cores used. Default doesn't use parallelism.
"_PACKAGE"
