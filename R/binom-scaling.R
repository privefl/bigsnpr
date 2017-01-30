################################################################################

#' Binomial(n, p) scaling
#'
#' Binomial(n, p) scaling where `n` is fixed and `p` is estimated.
#'
#' @param n Number of trials, parameter of the binomial distribution.
#' Default is `2`, which corresponds to diploidy, such as for the human genome.
#'
#' @return A new __function__ that returns a data.frame of two vectors
#' "mean" and "sd" which are of the length of __`ind.col`__.
#'
#' @details You will probably not use this function as is
#' but as the `fun.scaling` parameter of other functions of package `bigstatsr`.
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
#' X.svd <- bigstatsr::big_SVD(X, fun.scaling = snp_scaleBinom())
#' str(X.svd)
#' plot(X.svd$means)
#' abline(h = 2 * p, col = "red")
#' plot(X.svd$sds)
#' abline(h = sqrt(2 * p * (1 - p)), col = "red")
snp_scaleBinom <- function(n = 2) {
  function(X, ind.train = seq(nrow(X)), ind.col = seq(ncol(X))) {
    means <- bigstatsr::big_colstats(X, ind.train, ind.col)$sum /
      length(ind.train)
    p <- means / 2
    sds <- sqrt(2 * p * (1 - p))
    data.frame(mean = means, sd = sds)
  }
}

################################################################################

#' MAF
#'
#' Minor Allele Frequency.
#'
#' @param X The slot "genotypes" of a "bigSNP",
#' a [big.matrix][bigmemory::big.matrix-class].
#' You shouldn't have missing values in your data.
#'
#'
#' @param ind.train An optional vector of the row (individuals) indices that
#' are used, for the training part. If not specified, all rows are used.
#'
#' @param ind.col An optional vector of the column (SNP) indices that are used.
#' If not specified, all columns are used.
#'
#' @inheritParams snp_scaleBinom
#'
#' @return A vector of MAFs, correspond to `ind.col`.
#' @export
#'
snp_MAF <- function(X, ind.train = seq(nrow(X)),
                    ind.col = seq(ncol(X)), n = 2) {
  p <- bigstatsr::big_colstats(X, ind.train, ind.col)$sum /
    (n * length(ind.train))
  pmin(p, 1 - p)
}
