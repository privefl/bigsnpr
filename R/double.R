################################################################################

#' Double a genotype matrix
#'
#' Recode each column of a genotype matrix to be two columns, one encoding a
#' recessive status and the other a dominant status.
#'
#' @inheritParams bigsnpr-package
#'
#' @return The new genotype matrix, as a `big.matrix`.
#' @export
#'
#' @examples
snp_double <- function(x) {
  check_x(x)
  X <- x$genotypes

  newfile <- checkFile(x, "double")
  X2 <- big.matrix(nrow(X), 2 * ncol(X), type = "char",
                   backingfile = paste0(newfile, ".bk"),
                   backingpath = x$backingpath,
                   descriptorfile = paste0(newfile, ".desc"))

  doubleBM(X@address, X2@address)

  X2
}

################################################################################
