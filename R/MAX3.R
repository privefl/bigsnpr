################################################################################

ZCATT <- function(counts, val) {
  rj <- counts$cols.cases
  sj <- counts$cols.controls
  r <- colSums(rj)
  s <- colSums(sj)
  n <- r + s
  phi <- r / n
  p <- sweep(counts$cols.cases, 2, n, '/')
  q <- sweep(counts$cols.controls, 2, n, '/')

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
#'
#' @param val Computing `ZCATT(x)` for `x` in `val`. Default is `c(0, 0.5, 1)`
#' and corresponds to the _MAX3_ statistic. Only `c(0, 1)` corresponds to
#' _MAX2_. And only `0.5` corresponds to the Armatige trend test.
#' Finally, `seq(0, 1, length.out = L)` corresponds to _MAXL_.
#'
#' @examples
#' set.seed(1)
#'
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
#' fake$fam$affection <- c(rep(1, 5), rep(-1, 5))
#'
#' # Get MAX3 statistics
#' print(snp_MAX3(fake))
#'
#' @return A named list of __`S`__ and __`pS`__ for every column,
#' which are MAX3 statistics and associated p-values. __P-values are in
#' fact the minimum of the 3 p-values of each test separately.__ One can use
#' genomic control to rescale these p-values.
#'
#' @references Zheng, G., Yang, Y., Zhu, X., & Elston, R. (2012).
#' Robust Procedures. Analysis Of Genetic Association Studies, 151-206.
#' \url{http://dx.doi.org/10.1007/978-1-4614-2245-7_6}.
#'
#' @export
snp_MAX3 <- function(x, ind.train = seq(nrow(X)), val = c(0, 0.5, 1)) {
  check_x(x)

  X <- x$genotypes
  y <- transform_levels(x$fam$affection)

  cases <- (y == 1)

  ind.cases <- intersect(ind.train, which(cases))
  ind.controls <- intersect(ind.train, which(!cases))

  counts <- mycount2(X@address, ind.cases, ind.controls)

  stats <- ZCATT(counts, val)
  stats <- replace(stats, is.na(stats), 0)
  S <- apply(stats^2, 1, max)

  list(S = S, pS = pchisq(S, df = 1, lower.tail = FALSE))
}

################################################################################
