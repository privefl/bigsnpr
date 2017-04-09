################################################################################
#### Useful functions ####

rollMean <- function(x, size) {

  len <- 2 * size + 1

  lims <- qnorm(range(ppoints(len)))
  weights <- dnorm(seq(lims[1], lims[2], length.out = len))

  roll_mean(x, weights)
}

getIntervals <- function(x, n = 2) {

  le <- rle(diff(x))

  ind <- which((le$values == 1) & (le$lengths >= (n - 1)))
  pos <- cumsum(le$lengths) + 1

  cbind(x[c(1, pos)[ind]], x[pos[ind]])
}

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

#' Title
#'
#' @param G
#' @param ind.row
#' @param ind.col
#' @param fun.scaling
#' @param thr.r2.init
#' @param thr.r2.step
#' @param size
#' @param k
#' @param roll.size
#' @param int.min.size
#' @param ncores
#'
#' @return
#' @export
#'
#' @import foreach
#'
#' @examples
snp_clumpedSVD <- function(G,
                           ind.row = rows_along(G),
                           ind.col = cols_along(G),
                           fun.scaling = snp_scaleBinom(),
                           thr.r2.init = 0.2,
                           thr.r2.step = 2,
                           size = 500,
                           k = 10,
                           roll.size = 50,
                           int.min.size = 20,
                           ncores = 1) {

  # get BM
  G2 <- attach.BM(G)

  # first clumping
  THR <- thr.r2.init
  ind.keep <- snp_clumping(G, popres$map$chromosome,
                           thr.r2 = THR, size = size, ncores = ncores)

  repeat {
    # SVD
    obj.svd <- big_randomSVD(G,
                             fun.scaling = fun.scaling,
                             ncores = ncores,
                             ind.col = ind.keep,
                             k = k)

    # -log p-values of begin an outlier (by PC)
    lpval <- -2 * apply(abs(obj.svd$v), 2, function(x) {
      stats::pnorm(x, sd = mad(x), lower.tail = FALSE, log.p = TRUE)
    })
    # threshold of being an outlier based on Tukey's rule
    # http://math.stackexchange.com/a/966337
    lim <- stats::quantile(lpval, 0.75) + 1.5 * stats::IQR(lpval)

    # roll mean to get only consecutive outliers
    lpval.roll <- apply(lpval, 2, rollMean, size = roll.size)
    # get indices by intervals
    ind <- sort(unique(which(lpval.roll > lim, arr.ind = TRUE)[, "row"]))
    ind.range <- getIntervals(ind, int.min.size)

    if (nrow(ind.range)) {
      # local clumping on previous intervals
      THR <- THR / thr.r2.step
      ind.excl <- foreach(ic = rows_along(ind.range)) %do% {
        clumping.local(G2, ind.row, ind.keep[seq2(ind.range[ic, ])], THR)
      }
      ind.keep <- setdiff(ind.keep, unlist(ind.excl))
    } else { # converged
      break
    }
  }

  obj.svd
}

################################################################################
