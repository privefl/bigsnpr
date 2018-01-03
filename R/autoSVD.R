################################################################################
#### Useful functions ####

# apply a gaussian smoothing
rollMean <- function(x, size) {

  len <- 2 * size + 1

  lims <- stats::qnorm(range(stats::ppoints(len)))
  weights <- stats::dnorm(seq(lims[1], lims[2], length.out = len))

  roll_mean(x, weights)
}

# regroup consecutive integers in intervals
getIntervals <- function(x, n = 2) {

  le <- rle(diff(x))

  ind <- which((le$values == 1) & (le$lengths >= (n - 1)))
  pos <- cumsum(le$lengths) + 1

  cbind(x[c(1, pos)[ind]], x[pos[ind]])
}

################################################################################

#' Truncated SVD with pruning
#'
#' Fast truncated SVD which iteratively try to remove long-range LD regions
#' which appear in loadings of SVD.
#'
#' Improvements will come in the future, as for example, warm starts in order
#' to make the SVD computations faster.
#'
#' @inheritParams bigsnpr-package
#' @param fun.scaling A function that returns a named list of
#' `mean` and `sd` for every column, to scale each of their elements
#' such as followed: \deqn{\frac{X_{i,j} - mean_j}{sd_j}.}
#' Default is `snp_scaleBinom()`.
#' @param thr.r2 Threshold on squared correlation between SNPs.
#' Default is `0.2`.
#' @param size Radius of the window's size for the LD evaluations of the initial
#' step of clumping. Default is `500`.
#' @param k Number of singular vectors/values to compute.
#' Default is `10`. **This algorithm should be used to compute only
#' a few singular vectors/values.**
#' @param roll.size Radius of rolling windows to smooth log-p-values.
#'   Default is `50.`
#' @param int.min.size Minimum size of intervals of consecutive significant
#' indices to be reported as long-range LD region. Default is `20`.
#' @param verbose Output some information on the iterations? Default is `TRUE`.
#'
#' @inherit bigstatsr::big_randomSVD return
#' @export
#'
#' @import foreach
#'
#' @examples
#' ex <- snp_attachExtdata()
#'
#' obj.svd <- snp_autoSVD(G = ex$genotypes,
#'                        infos.chr = ex$map$chromosome,
#'                        infos.pos = ex$map$physical.position)
#'
#' str(obj.svd)
#'
snp_autoSVD <- function(G,
                        infos.chr,
                        infos.pos,
                        ind.row = rows_along(G),
                        ind.col = cols_along(G),
                        fun.scaling = snp_scaleBinom(),
                        thr.r2 = 0.2,
                        size = 500,
                        k = 10,
                        roll.size = 50,
                        int.min.size = 20,
                        ncores = 1,
                        verbose = TRUE) {

  check_args()

  # Verbose?
  printf2 <- function(...) if (verbose) printf(...)

  # First clumping
  THR <- thr.r2
  printf2("Phase of clumping at r2 > %s.. ", THR)
  ind.keep <- snp_clumping(G, infos.chr,
                           thr.r2 = THR,
                           size = size,
                           ncores = ncores)
  printf2("keep %d SNPs.\n", length(ind.keep))

  iter <- 1
  LRLDR <- LD.wiki34[0, 1:3]
  repeat {
    printf2("\nIteration %d:\n", iter)
    printf2("Computing SVD..\n")
    # SVD
    obj.svd <- big_randomSVD(G,
                             fun.scaling = fun.scaling,
                             ncores = ncores,
                             ind.col = ind.keep,
                             k = k)

    # The -log10 p-values of being an outlier (by PC)
    lpval <- -stats::predict(snp_pcadapt(G, obj.svd$u, ind.col = ind.keep))
    # Threshold of being an outlier based on Tukey's rule
    # http://math.stackexchange.com/a/966337
    lim <- stats::quantile(lpval, 0.75) + 1.5 * stats::IQR(lpval)

    # Roll mean to get only consecutive outliers and
    ind.excl <- which(rollMean(lpval, size = roll.size) > lim)
    # Regroup them by intervals to return long-range LD regions
    ind.range <- getIntervals(ind.excl, n = int.min.size)
    if (nrow(ind.range) > 0) {
      LRLDR.add <- cbind(
        infos.chr[ind.keep[ind.range[, 1]]],
        matrix(infos.pos[ind.keep[ind.range]], ncol = 2)
      )
      LRLDR[nrow(LRLDR) + rows_along(LRLDR.add), ] <- LRLDR.add
    }

    if (length(ind.excl) > 0) {
      printf2("%d long-range LD regions were detected..\n", nrow(ind.range))
      ind.keep <- ind.keep[-ind.excl]
      iter <- iter + 1
    } else { # converged
      printf2("\nConverged!\n")
      break
    }
  }

  structure(obj.svd, subset = ind.keep, lrldr = LRLDR)
}

################################################################################
