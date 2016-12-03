################################################################################

#' Counts
#'
#' Counts the number of 0, 1, 2
#' by SNP and phenotype (cases/controls) and the number
#' of NA by individual. The number of NA by SNP can be deduced
#' using the relation \eqn{n_{NA} = n - (n_0 + n_1 + n_2)}.
#'
#' @inheritParams bigsnpr-package
#' @examples
#' set.seed(1)
#'
#' # constructing a fake genotype big.matrix
#' a <- big.matrix(10, 15, type = "char")
#' a[] <- sample(c(0, 1, 2, NA), 150, TRUE)
#' print(a[,])
#'
#' # constructing a fake incomplete bigSNP with 10 individuals and 15 SNPs
#' # where the 5 first individuals are cases and the 5 last are controls.
#' fake <- list()
#' class(fake) <- "bigSNP"
#' fake$genotypes <- a
#' fake$fam$affection <- c(rep(2, 5), rep(1, 5))
#'
#' # Get counts
#' print(test1 <- snp_counts(fake))
#' print(test2 <- snp_counts(fake, has.pheno = FALSE))
#'
#' # Same results
#' test1$cols <- test1$cols.controls + test1$cols.cases
#' print(all.equal(test1$cols, test2$cols))
#' print(all.equal(test1$rows, test2$rows))
#'
#' @return A list of:
#' - one or two matrices of size 3*m (m is the number of SNPs)
#' representing the counts of 0, 1 and 2 (by status),
#' - an integer vector of size n (n is the number of individuals)
#' corresponding to the number of missing values per individual.
#' @export
snp_counts <- function(x, has.pheno = TRUE) {
  X <- x$genotypes
  if (has.pheno) {
    y <- transform_levels(x$fam$affection)
    ind.cases    <- which(y == 1)
    ind.controls <- which(y == 0)
    res <- mycount2(X@address, ind.cases, ind.controls)
  } else {
    res <- mycount2(X@address, integer(0), seq(nrow(X)))
    res$cols.cases <- NULL
    names(res)[1] <- "cols"
  }

  res
}

################################################################################
