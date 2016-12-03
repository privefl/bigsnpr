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

#' MAX3 statistic
#'
#' Counts the number of 0, 1, 2 by SNP and phenotype (cases/controls)
#' and then compute the MAX3 statistic.
#'
#' @inheritParams bigsnpr-package
#' @examples
#' set.seed(1)
#'
#' # constructing a fake genotype big.matrix
#' a <- big.matrix(10, 15, type = "char")
#' a[] <- sample(c(0, 1, 2), 150, TRUE)
#' print(a[,])
#'
#' # constructing a fake incomplete bigSNP with 10 individuals and 15 SNPs
#' # where the 5 first individuals are cases and the 5 last are controls.
#' fake <- list()
#' class(fake) <- "bigSNP"
#' fake$genotypes <- a
#' fake$fam$affection <- c(rep(1, 5), rep(-1, 5))
#'
#' # Get MAX3 statistics
#' print(snp_MAX3(fake))
#'
#' @return A named list of __`S`__ and __`pS`__ for every column,
#' which are MAX3 statistics and associated p-values. __P-values are in
#' fact the minimum of the 3 p-values of each test separately.__ One can use
#' genomic control to rescale these p-values.
#' @export
snp_MAX3 <- function(x, ind.train = seq(nrow(X))) {
  check_x(x)

  X <- x$genotypes
  y <- transform_levels(x$fam$affection)

  cases <- (y == 1)

  ind.cases <- intersect(ind.train, which(cases))
  ind.controls <- intersect(ind.train, which(!cases))

  counts <- mycount(X@address, ind.cases, ind.controls)

  stats <- sapply(c(0, 0.5, 1), function(x) ZCATT(counts, x = x))
  stats <- replace(stats, is.na(stats), 0)
  S <- apply(stats^2, 1, max)

  list(S = S, pS = pchisq(S, df = 1, lower.tail = FALSE))
}

################################################################################
