################################################################################

#' @import bigstatsr
#' @useDynLib bigsnpr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#'
#' @param G A [FBM.code256][FBM.code256-class]
#' (typically `<bigSNP>$genotypes`).\cr
#' **You shouldn't have missing values in your data or SNPs with 0 MAF.**
#' @param Gna A [FBM.code256][FBM.code256-class]
#' (typically `<bigSNP>$genotypes`).\cr
#' You can have missing values in your data.
#'
#' @param x A [bigSNP][bigSNP-class].
#'
#' @param infos.chr Vector of integers specifying each SNP's chromosome.\cr
#' Typically `<bigSNP>$map$chromosome`.
#'
#' @param infos.pos Vector of integers specifying the physical position
#' on a chromosome (in base pairs) of each SNP.\cr
#' Typically `<bigSNP>$map$physical.pos`.
#'
#' @param nploidy Number of trials, parameter of the binomial distribution.
#' Default is `2`, which corresponds to diploidy, such as for the human genome.
#'
#' @param ind.row An optional vector of the row indices (individuals) that
#' are used. If not specified, all rows are used.\cr
#' **Don't use negative indices.**
#'
#' @param ind.col An optional vector of the column indices (SNPs) that are used.
#' If not specified, all columns are used.\cr
#' **Don't use negative indices.**
#'
#' @param ncores Number of cores used. Default doesn't use parallelism.
#'   You may use [nb_cores].
#'
"_PACKAGE"

################################################################################
