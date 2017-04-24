################################################################################

ZCATT <- function(counts.cases, counts.controls, val) {
  rj <- counts.cases
  sj <- counts.controls
  r <- colSums(rj)
  s <- colSums(sj)
  n <- r + s
  phi <- r / n
  p <- sweep(counts.cases, 2, n, '/')
  q <- sweep(counts.controls, 2, n, '/')

  num <- sweep(rj, 2, 1 - phi, '*') - sweep(sj, 2, phi, '*')
  pj <- sweep(rj + sj, 2, n, '/')
  coef <- n * phi * (1 - phi)

  ZCATT.part <- function(x) {
    x2 <- c(0, x, 1)

    num2 <- colSums(x2 * num)
    deno <- colSums(x2^2*pj) - (colSums(x2*pj))^2
    deno2 <- sqrt(coef * deno)

    num2 / deno2
  }

  sapply(val, ZCATT.part)
}

################################################################################

#' MAX3 statistic
#'
#' Counts the number of 0, 1, 2 by SNP and phenotype (cases/controls)
#' and then compute the MAX3 statistic.
#'
#' @inheritParams bigsnpr-package
#' @param ind.train An optional vector of the row indices that are used,
#' for the training part. If not specified, all rows are used.\cr
#' __Don't use negative indices.__
#' @param val
#' Computing \eqn{\smash{\displaystyle\max_{x \in val}}~Z_{CATT}^2(x)}.
#' Default is `c(0, 0.5, 1)` and corresponds to the _MAX3_ statistic.
#' Only `c(0, 1)` corresponds to _MAX2_.
#' And only `0.5` corresponds to the Armitage trend test.
#' Finally, `seq(0, 1, length.out = L)` corresponds to _MAXL_.
#'
#' @examples
#' set.seed(1)
#'
#' # constructing a fake genotype big.matrix
#' a <- big.matrix(10, 12, type = "char", shared = FALSE)
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
#' @return A data.frame of `S` and `pS` for every column, which are MAX3
#' statistics and associated p-values. __P-values are in fact the minimum of
#' the p-values of each test separately (so they are biased downward).__
#' One can use genomic control to rescale these p-values.
#'
#' @references Zheng, G., Yang, Y., Zhu, X., & Elston, R. (2012).
#' Robust Procedures. Analysis Of Genetic Association Studies, 151-206.
#' \url{http://dx.doi.org/10.1007/978-1-4614-2245-7_6}.
#'
#' @export
snp_MAX3 <- function(x, ind.train = rows_along(X.), val = c(0, 0.5, 1)) {

  X. <- x$genotypes
  y <- transform_levels(x$fam$affection)

  cases <- (y == 1)

  ind.cases <- intersect(ind.train, which(cases))
  ind.controls <- intersect(ind.train, which(!cases))

  onlyO12 <- as.character(0:2) # don't consider missing values
  counts.cases <- big_counts(X., ind.row = ind.cases)[onlyO12, ]
  counts.controls <- big_counts(X., ind.row = ind.controls)[onlyO12, ]

  stats <- ZCATT(counts.cases, counts.controls, val)
  stats <- replace(stats, is.na(stats), 0)

  S <- apply(stats^2, 1, max)

  data.frame(S = S, pS = stats::pchisq(S, df = 1, lower.tail = FALSE))
}

################################################################################
