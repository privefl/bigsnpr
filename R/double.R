################################################################################

#' Double a genotype matrix
#'
#' Recode each column of a genotype matrix to be two columns, one encoding a
#' recessive status and the other a dominant status. Values will be rounded
#' so that you have only 0s, 1s and 2s in your data.
#'
#' @inheritParams bigsnpr-package
#'
#' @return The new "doubled" genotype matrix, as a `BM.code.descriptor`,
#' which values are 0s and 1s.
#'
#' @examples
#' test <- snp_attachExtdata()
#' X <- attach.BM(test$genotypes)
#' X[1:8, 1:5]
#'
#' test2 <- snp_double(test)
#' attach.BM(test2)[1:8, 1:10]
#'
#' rm(X)
#'
#' @export
snp_double <- function(x) {

  X <- attach.BM(x$genotypes)
  X@code <- round(CODE_DOSAGE)

  newfiles <- getNewFiles(x$savedIn, "doubled")

  X2 <- big.matrix(nrow(X), 2 * ncol(X), type = "raw",
                   backingfile = basename(newfiles$bk),
                   backingpath = dirname(newfiles$bk),
                   descriptorfile = basename(newfiles$desc))

  doubleBM(X, X2@address)

  describe(as.BM.code(X2, code = CODE_01))
}

################################################################################
