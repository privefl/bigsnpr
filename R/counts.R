################################################################################

#' @title Counts.
#' @name counts
#' @inheritParams bigsnpr-package
#' @examples
#' # constructing a fake genotype big.matrix
#' a <- big.matrix(10, 15, type = "char", shared = FALSE)
#' a[] <- sample(c(0, 1, 2, NA), 150, TRUE)
#' print(a[,])
#'
#' # constructing a fake incomplete bigSNP with 10 individuals and 15 SNPs
#' # where the 5 first individuals are cases and the 5 last are controls.
#' fake <- list()
#' class(fake) <- "bigSNP"
#' fake$genotypes <- a
#' fake$fam$pheno <- c(rep(1, 5), rep(-1, 5))
#'
#' # Get counts
#' print(CountByPheno(fake))
#' print(CountNAByRow(fake))
NULL

################################################################################

#' @description \code{CountByPheno}: counts the number of 0, 1, 2 and NA
#' by SNP and phenotype (cases/controls).
#' @return \code{CountByPheno}: An integer matrix
#' of size 8*m (m is the number of SNPs).
#' @rdname counts
#' @export
CountByPheno <- function(x) {
  ind.cases <- which(x$fam$pheno == 1)
  ind.controls <- which(x$fam$pheno == -1)
  res <- mycount2((x$genotypes)@address, ind.cases, ind.controls)
  row.names(res) <- paste0("nb.", c(paste0(c(0:2, NA), ".case"),
                                    paste0(c(0:2, NA), ".control")))

  return(res)
}

################################################################################

#' @description \code{CountNAByRow}: counts the number of NA
#' by individuals.
#' @return \code{CountNAByRow}: An integer vector
#' of size n (n is the number of individuals).
#' @rdname counts
#' @export
CountNAByRow <- function(x) {
  mycount3((x$genotypes)@address)
}

################################################################################

