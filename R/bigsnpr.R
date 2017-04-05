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
#' @param has.pheno Does __`x`__ has well-specified phenotypes?
#' i.e. `x$fam$affection` is a vector with exactly two values
#' (and the value coding for cases is larger than the value coding
#' for controls). So, in particular, no missing value is allowed.
#' Default is `TRUE`.
#'
#' @param nploidy Number of trials, parameter of the binomial distribution.
#' Default is `2`, which corresponds to diploidy, such as for the human genome.
#'
#' @param ind.row An optional vector of the row (individuals) indices that
#' are used. If not specified, all rows are used.
#'
#' @param ind.col An optional vector of the column (SNP) indices that are used.
#' If not specified, all columns are used.
#'
#' @param block.size Maximum number of loci read at once (for all individuals).
#'
#' @param ncores Number of cores used. Default doesn't use parallelism.
"_PACKAGE"

################################################################################
