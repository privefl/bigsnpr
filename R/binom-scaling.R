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
snp_scaleBinom <- function(nploidy = 2) {

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
                    nploidy = 2) {

  p <- big_colstats(G, ind.row = ind.row, ind.col = ind.col)$sum /
    (nploidy * length(ind.row))

  pmin(p, 1 - p)
}

################################################################################

#' Binomial(2, p) scaling
#'
#' Binomial(2, p) scaling where `p` is estimated.
#'
#' @inheritParams bed_autoSVD
#'
#' @return A data frame with `$center` and `$scale`.
#'
#' @details You will probably not use this function as is but as the
#' `fun.scaling` parameter of other functions of package `bigstatsr`.
#'
#' @export
#'
#' @inherit snp_scaleBinom references
#'
#' @examples
#' bedfile <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
#' obj.bed <- bed(bedfile)
#'
#' str(bed_scaleBinom(obj.bed))
#'
#' str(bed_randomSVD(obj.bed, bed_scaleBinom))
#'
bed_scaleBinom <- function(obj.bed,
                           ind.row = rows_along(obj.bed),
                           ind.col = cols_along(obj.bed)) {

  stats <- bed_stats(obj.bed, ind.row, ind.col)
  af <- stats$sum / (2 * stats$nb_nona_col)

  data.frame(center = 2 * af, scale = sqrt(2 * af * (1 - af)))
}

################################################################################

#' Counts
#'
#' Counts the number of 0s, 1s, 2s and NAs by variants in the bed file.
#'
#' @inheritParams bigsnpr-package
#'
#' @return A matrix of with 4 rows and `length(ind.col)` columns.
#'
#' @export
#'
#' @examples
#' bedfile <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
#' obj.bed <- bed(bedfile)
#'
#' bed_counts(obj.bed, ind.col = 1:5)
#'
bed_counts <- function(obj.bed,
                       ind.row = rows_along(obj.bed),
                       ind.col = cols_along(obj.bed),
                       ncores = 1) {

  res <- big_parallelize(obj.bed, p.FUN = function(X, ind, ind.row) {
    bed_counts_cpp(obj.bed, ind.row, ind)
  }, p.combine = "cbind", ncores = ncores, ind = ind.col, ind.row = ind.row)

  rownames(res) <- c(0:2, NA)
  res
}

################################################################################

#' Allele frequencies
#'
#' Allele frequencies of a [bed] object.
#'
#' @inheritParams bigsnpr-package
#'
#' @return A data.frame with
#'  - `$ac`: allele counts,
#'  - `$mac`: minor allele counts,
#'  - `$af`: allele frequencies,
#'  - `$maf`: minor allele frequencies,
#'  - `$N`: numbers of non-missing values.
#'
#' @export
#'
#' @examples
#' bedfile <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
#' obj.bed <- bed(bedfile)
#'
#' bed_MAF(obj.bed, ind.col = 1:5)
#'
bed_MAF <- function(obj.bed,
                    ind.row = rows_along(obj.bed),
                    ind.col = cols_along(obj.bed),
                    ncores = 1) {

  counts <- bed_counts(obj.bed, ind.row, ind.col, ncores)
  ac <- counts[2, ] + 2 * counts[3, ]
  nb_nona <- length(ind.row) - counts[4, ]
  af <- ac / (2 * nb_nona)

  data.frame(
    ac  = ac,
    mac = pmin(ac, 2 * nb_nona - ac),
    af  = af,
    maf = pmin(af, 1 - af),
    N   = nb_nona
  )
}

################################################################################


