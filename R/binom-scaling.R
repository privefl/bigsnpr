################################################################################

#' Binomial(n, p) scaling
#'
#' Binomial(n, p) scaling where `n` is fixed and `p` is estimated.
#'
#' @inheritParams bigsnpr-package
#'
#' @return A new __function__ that returns a data.frame of two vectors
#' "center" and "scale" which are of the length of `ind.col`.
#'
#' @details You will probably not use this function as is but as the
#' `fun.scaling` parameter of other functions of package `bigstatsr`.
#'
#' @export
#'
#' @references This scaling is widely used for SNP arrays.
#' Patterson N, Price AL, Reich D (2006).
#' Population Structure and Eigenanalysis.
#' PLoS Genet 2(12): e190.
#' \url{http://dx.doi.org/10.1371/journal.pgen.0020190}.
#'
#' @examples
#' set.seed(1)
#'
#' a <- matrix(0, 93, 170)
#' p <- 0.2
#' a[] <- rbinom(length(a), 2, p)
#' X <- add_code256(big_copy(a, type = "raw"), code = c(0, 1, 2, rep(NA, 253)))
#' X.svd <- big_SVD(X, fun.scaling = snp_scaleBinom())
#' str(X.svd)
#' plot(X.svd$center)
#' abline(h = 2 * p, col = "red")
#' plot(X.svd$scale)
#' abline(h = sqrt(2 * p * (1 - p)), col = "red")
snp_scaleBinom <- function(nploidy = getOption("bigsnpr.nploidy")) {

  function(X,
           ind.row = rows_along(X),
           ind.col = cols_along(X)) {

    means <- big_colstats(X, ind.row = ind.row, ind.col = ind.col)$sum /
      length(ind.row)

    p <- means / nploidy
    sds <- sqrt(nploidy * p * (1 - p))

    data.frame(center = means, scale = sds)
  }
}

################################################################################

#' MAF
#'
#' Minor Allele Frequency.
#'
#' @inheritParams bigsnpr-package
#'
#' @return A vector of MAFs, corresponding to `ind.col`.
#' @export
#'
#' @examples
#' obj.bigsnp <- snp_attachExtdata()
#' str(maf <- snp_MAF(obj.bigsnp$genotypes))
#'
snp_MAF <- function(G,
                    ind.row = rows_along(G),
                    ind.col = cols_along(G),
                    nploidy = getOption("bigsnpr.nploidy")) {

  p <- big_colstats(G, ind.row = ind.row, ind.col = ind.col)$sum /
    (nploidy * length(ind.row))

  pmin(p, 1 - p)
}

################################################################################
