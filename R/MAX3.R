################################################################################

ZCATT <- function(counts, x) {
  r <- sum(counts$cases[, 1])
  p <- counts$cases / r
  s <- sum(counts$controls[, 1])
  q <- counts$controls / s
  x <- c(0, x, 1)
  num <- colSums(x*(p - q))
  deno1 <- colSums(x^2*p) - (colSums(x*p))^2
  deno2 <- colSums(x^2*q) - (colSums(x*q))^2
  deno <- sqrt(deno1/r + deno2/s)
  num / deno
}

################################################################################

#' @title MAX3 statistic.
#' @name MAX3
#' @inheritParams bigsnpr-package
#' @examples
#' # constructing a fake genotype big.matrix
#' a <- big.matrix(10, 15, type = "char", shared = FALSE)
#' a[] <- sample(c(0, 1, 2), 150, TRUE)
#' print(a[,])
#'
#' # constructing a fake incomplete bigSNP with 10 individuals and 15 SNPs
#' # where the 5 first individuals are cases and the 5 last are controls.
#' fake <- list()
#' class(fake) <- "bigSNP"
#' fake$genotypes <- a
#' fake$fam$pheno <- c(rep(1, 5), rep(-1, 5))
#'
#' # Get MAX3 statistics
#' print(MAX3(fake))
#' @description Counts the number of 0, 1, 2
#' by SNP and phenotype (cases/controls) and then
#' compute the MAX3 statistic.
#' @return A numeric vector of the column statistics.
#' @export
MAX3 <- function(x, ind.train = seq(nrow(x$genotypes))) {
  check_x(x)

  cases <- (x$fam$pheno == 1)

  ind.cases <- intersect(ind.train, which(cases))
  ind.controls <- intersect(ind.train, which(!cases))

  counts <- mycount((x$genotypes)@address, ind.cases, ind.controls)

  stats <- sapply(c(0, 0.5, 1), function(x) ZCATT(counts, x = x))
  stats <- replace(stats, is.na(stats), 0)

  apply(abs(stats), 1, max)
}

################################################################################
