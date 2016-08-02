#' @import bigmemory
#' @useDynLib bigsnpr
#' @importFrom Rcpp sourceCpp
#' @title R Package for analysis of massive SNP arrays.
#' @description TODO
#' @name bigsnpr-package
#' @param x A \code{bigSNP}.
#' @param ind.train An optional vector of the row indices that are used,
#' for the training part.
#' If not specified, all data are used.
#' @param ncores Number or cores used.
#' Default doesn't use parallelism.
#' @aliases bigsnpr-package bigsnpr
#' @keywords package
NULL
