################################################################################

#' @title Counts.
#' @name Counts
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
#' print(Counts(fake))
#' @description Counts the number of 0, 1, 2
#' by SNP and phenotype (cases/controls) and the number
#' of NA by individual. The number of NA by SNP can be deduced
#' using the relation \eqn{n_{NA} = n - (n_0 + n_1 + n_2)}.
#' @return A list of:\itemize{
#' \item an integer matrix of size 6*m (m is the number of SNPs),
#' \item an integer vector of size n (n is the number of individuals).
#' }
#' @export
Counts <- function(x) {
  ind.cases <- which(x$fam$pheno == 1)
  ind.controls <- which(x$fam$pheno == -1)
  res <- mycount2((x$genotypes)@address, ind.cases, ind.controls)

  return(res)
}

################################################################################
