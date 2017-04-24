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

# get indices to exclude in a small region of LD
clumping.local <- function(G2, ind.row, ind.col, thr.r2) {

  # cache some computations
  stats <- big_colstats(G2, ind.row = ind.row, ind.col = ind.col)
  n <- length(ind.row)
  denoX <- (n - 1) * stats$var
  nploidy <- getOption("bigsnpr.nploidy")
  p <- stats$sum / (nploidy * n)
  maf <- pmin(p, 1 - p)

  # main algo
  keep <- local_clumping(G2,
                         rowInd = ind.row,
                         colInd = ind.col,
                         ordInd = order(maf, decreasing = TRUE),
                         sumX = stats$sum,
                         denoX = denoX,
                         thr = thr.r2)

  ind.col[!keep]
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
#' @param thr.r2.init Threshold on squared correlation between SNPs.
#' Default is `0.2`.
#' @param thr.r2.step Coefficient of division of `thr.r2` at each step.
#' Default is `2` which means that at iteration 2, `thr.r2 = 0.1`, at iteration
#' 3, `thr.r2 = 0.05`, etc.
#' @param size Radius of the window's size for the LD evaluations of the initial
#' step of clumping. Default is `500`.
#' @param k Number of singular vectors/values to compute.
#' Default is `10`. **This algorithm should be used to compute only
#' a few singular vectors/values.**
#' @param roll.size Radius of rolling windows to smooth log-p-values.
#' @param int.min.size Minimum size of intervals of consecutive significant
#' indices.
#' @param verbose Output some information on the iterations? Default is `TRUE`.
#'
#' @inherit bigstatsr::big_randomSVD return
#' @export
#'
#' @import foreach
#' @importFrom magrittr %>%
#'
#' @examples
#' ex <- snp_attachExtdata()
#'
#' obj.svd <- snp_clumpedSVD(G = ex$genotypes,
#'                           infos.chr = ex$map$chromosome)
#'
#' str(obj.svd)
#'
snp_clumpedSVD <- function(G,
                           infos.chr,
                           ind.row = rows_along(G),
                           ind.col = cols_along(G),
                           fun.scaling = snp_scaleBinom(),
                           thr.r2.init = 0.2,
                           thr.r2.step = 2,
                           size = 500,
                           k = 10,
                           roll.size = 50,
                           int.min.size = 20,
                           ncores = 1,
                           verbose = TRUE) {

  # get BM
  G2 <- attach.BM(G)

  # verbose?
  printf2 <- function(...) if (verbose) printf(...)

  # first clumping
  THR <- thr.r2.init
  printf2("First phase of clumping at r2 > %s.. ", THR)
  ind.keep <- snp_clumping(G, infos.chr,
                           thr.r2 = THR,
                           size = size,
                           ncores = ncores)
  printf2("keep %d SNPs.\n", length(ind.keep))

  iter <- 1
  repeat {
    printf2("\nIteration %d:\n", iter)
    printf2("Computing SVD..\n")
    # SVD
    obj.svd <- big_randomSVD(G,
                             fun.scaling = fun.scaling,
                             ncores = ncores,
                             ind.col = ind.keep,
                             k = k)

    # -log p-values of being an outlier (by PC)
    lpval <- -2 * apply(abs(obj.svd$v), 2, function(x) {
      stats::pnorm(x, sd = stats::mad(x), lower.tail = FALSE, log.p = TRUE)
    })
    # threshold of being an outlier based on Tukey's rule
    # http://math.stackexchange.com/a/966337
    lim <- stats::quantile(lpval, 0.75) + 1.5 * stats::IQR(lpval)

    # roll mean to get only consecutive outliers and regroup them by intervals
    ind.range <-
      apply(lpval, 2, rollMean, size = roll.size) %>%
      which(. > lim, arr.ind = TRUE)[, "row"] %>%
      unique() %>%
      sort() %>%
      getIntervals(n = int.min.size)


    if (N <- nrow(ind.range)) {
      # local clumping on previous intervals
      THR <- THR / thr.r2.step
      printf2("Local clumping in %d regions at r2 > %s.. ", N, THR)
      ind.excl <- foreach(ic = rows_along(ind.range), .combine = 'c') %do% {
        clumping.local(G2, ind.row, ind.keep[seq2(ind.range[ic, ])], THR)
      }
      ind.keep <- setdiff(ind.keep, ind.excl)
      printf2("further excluded %d SNPs.\n", length(ind.excl))
      iter <- iter + 1
    } else { # converged
      printf2("\nConverged!\n")
      break
    }
  }

  structure(obj.svd, subset = ind.keep)
}

################################################################################
