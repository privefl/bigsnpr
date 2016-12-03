#' Binomial(2, p) scaling
#'
#' @inheritParams bigsnpr-package
#'
#' @return A named list of two vectors __`mean`__ and __`sd`__
#' which are as long as the number of columns of __`X`__.
#'
#' @details You will probably not use this function as is
#' but as the __`fun.scaling`__ parameter of other functions
#' of package __`bigstatsr`__.
#' @export
#'
#' @references This scaling is widely used for SNP arrays.
#' Patterson N, Price AL, Reich D (2006)
#' Population Structure and Eigenanalysis.
#' PLoS Genet 2(12): e190.
#' \url{http://dx.doi.org/10.1371/journal.pgen.0020190}.
#'
#' @examples
#' set.seed(1)
#'
#' a <- matrix(0, 43, 17)
#' p <- 0.2
#' a[] <- rbinom(length(a), 2, p)
#' X <- as.big.matrix(a, type = "char", shared = FALSE)
#' X.svd <- bigstatsr::big_SVD(X, fun.scaling = snp_scaleBinom)
#' str(X.svd)
#' plot(X.svd$means)
#' abline(h = 2 * p, col = "red")
#' plot(X.svd$sds)
#' abline(h = sqrt(2 * p * (1 - p)), col = "red")
snp_scaleBinom <- function(X, ind.train = seq(nrow(X))) {
  means <- bigstatsr::big_colstats(X, ind.train)$sum /
    length(ind.train)
  p <- means / 2
  sds <- sqrt(2 * p * (1 - p))
  list(mean = means, sd = sds)
}
