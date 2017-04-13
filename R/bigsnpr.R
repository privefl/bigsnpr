################################################################################

#' @import bigmemory bigstatsr
#' @useDynLib bigsnpr
#' @importFrom Rcpp sourceCpp
#'
#' @param G A [BM.code.descriptor][BM.code.descriptor-class].
#' Typically the slot `genotypes` of a `bigSNP`.
#' You shouldn't have missing values in your data or SNPs with 0 MAF.
#'
#' @param x A [bigSNP][bigSNP-class].
#'
#' @param infos.chr Vector of integers specifying the chromosome of each SNP.
#' Typically the slot `chromosome` of the slot `map` of a `bigSNP`.
#'
#' @param infos.pos Vector of integers specifying the physical position
#' on a chromosome (in base pairs) of each SNP.
#' Typically the slot `physical.pos` of the slot `map` of a `bigSNP`.
#'
#' @param nploidy Number of trials, parameter of the binomial distribution.
#' Default is `2`, which corresponds to diploidy, such as for the human genome.
#'
#' @param ind.row An optional vector of the row indices (individuals) that
#' are used. If not specified, all rows are used.
#' **Don't use negative indices.**
#'
#' @param ind.col An optional vector of the column indices (SNPs) that are used.
#' If not specified, all columns are used.
#' **Don't use negative indices.**
#'
#' @param ncores Number of cores used. Default doesn't use parallelism.
"_PACKAGE"

################################################################################
