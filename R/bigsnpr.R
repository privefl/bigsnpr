#' @import bigmemory bigstatsr
#' @useDynLib bigsnpr
#' @importFrom Rcpp sourceCpp
#' @name bigsnpr-package
#' @param x A [bigSNP][bigSNP-class].
#' @param ind.train An optional vector of the row indices that are used,
#' for the training part. If not specified, all data are used.
#' @param block.size Maximum number of loci read at once (for all individuals).
#' @param ncores Number or cores used. Default doesn't use parallelism.
"_PACKAGE"
