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
#' Compute the MAX3 statistic, which tests for three genetic models
#' (additive, recessive and dominant).
#'
#' __P-values associated with returned scores are in fact the minimum of the
#' p-values of each test separately. Thus, they are biased downward.__
#'
#' @inheritParams bigsnpr-package
#' @inheritParams bigstatsr::`bigstatsr-package`
#' @param val
#' Computing \eqn{\smash{\displaystyle\max_{x \in val}}~Z_{CATT}^2(x)}.
#' - Default is `c(0, 0.5, 1)` and corresponds to the _MAX3_ statistic.
#' - Only `c(0, 1)` corresponds to _MAX2_.
#' - And only `0.5` corresponds to the Armitage trend test.
#' - Finally, `seq(0, 1, length.out = L)` corresponds to _MAXL_.
#'
#' @examples
#' set.seed(1)
#'
#' # constructing a fake genotype big.matrix
#' N <- 50; M <- 1200
#' fake <- snp_fake(N, M)
#' G <- fake$genotypes
#' G[] <- sample(as.raw(0:3), size = length(G), replace = TRUE)
#' G[1:8, 1:10]
#'
#' # Specify case/control phenotypes
#' fake$fam$affection <- rep(1:2, each = N / 2)
#'
#' # Get MAX3 statistics
#' y01 <- fake$fam$affection - 1
#' str(test <- snp_MAX3(fake$genotypes, y01.train = y01))
#' # p-values are not well calibrated
#' snp_qq(test)
#' # genomic control is not of much help
#' snp_qq(snp_gc(test))
#'
#' # Armitage trend test (well calibrated because only one test)
#' test2 <- snp_MAX3(fake$genotypes, y01.train = y01, val = 0.5)
#' snp_qq(test2)
#'
#' @inherit snp_pcadapt return
#'
#' @references Zheng, G., Yang, Y., Zhu, X., & Elston, R. (2012).
#' Robust Procedures. Analysis Of Genetic Association Studies, 151-206.
#' \url{http://dx.doi.org/10.1007/978-1-4614-2245-7_6}.
#'
#' @export
snp_MAX3 <- function(Gna, y01.train,
                     ind.train = rows_along(Gna),
                     val = c(0, 0.5, 1)) {

  check_args()

  is.case <- (y01.train == 1)
  ind.cases <- ind.train[is.case]
  ind.controls <- ind.train[!is.case]

  only012 <- as.character(0:2) # don't consider missing values
  counts.cases <- big_counts(Gna, ind.row = ind.cases)[only012, ]
  counts.controls <- big_counts(Gna, ind.row = ind.controls)[only012, ]

  stats <- ZCATT(counts.cases, counts.controls, val)
  stats <- replace(stats, is.na(stats), 0)

  fun.pred <- function(xtr) {
    stats::pchisq(xtr, df = 1, lower.tail = FALSE, log.p = TRUE) / log(10)
  }

  structure(data.frame(score = apply(stats^2, 1, max)),
            class = c("mhtest", "data.frame"),
            transfo = identity,
            predict = fun.pred)
}

################################################################################
