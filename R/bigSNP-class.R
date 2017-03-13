################################################################################

# code genotype calls or missing values
CODE_012 <- c(0, 1, 2, rep(NA, 253))

# used in `snp_double`
CODE_01  <- c(0, 1, rep(NA, 254))

# code genotype calls, missing values, imputed calls and imputed allele dosages
CODE_DOSAGE <- c(0, 1, 2, NA, 0, 1, 2, seq(0, 2, by = 0.01), rep(NA, 48))

################################################################################

#' S3 class bigSNP
#'
#' A class for representing infos on massive SNP arrays.
#'
#' A named list with at least 5 slots: \describe{
#'   \item{genotypes}{A filebacked [big.matrix][big.matrix-class]
#'     of type `char` (one byte signed integer), representing genotypes
#'     (read from a ".bed" file). Each element is either `0`, `1`, `2` or `NA`.
#'     Rows are individuals and columns are SNPs.}
#'   \item{fam}{A `data.frame` containing some information on the SNPs
#'     (read from a ".fam" file).}
#'   \item{map}{A `data.frame` giving some information on the individuals
#'     (read from a ".bim" file).}
#'   \item{backingfile}{The root name for the backing files of the object.}
#'   \item{backingpath}{The path to the directory containing the backing files.}
#' }
#'
#' @name bigSNP-class
#' @aliases bigSNP-class bigSNP
#' @keywords class
#' @seealso [snp_readBed]
NULL

################################################################################
