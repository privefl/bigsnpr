################################################################################
#### Useful functions ####

# Outlier detection (upper)
tukey_MC_up <- function(x, coef = NULL, a = -4, b = 3, alpha = 0.05) {

  if (is.null(coef)) {
    alpha <- alpha / sum(!is.na(x))  # Bonferroni correction
    coef <- (qnorm(alpha, lower.tail = FALSE) - 0.6744898) / 1.34898
  }

  robustbase::adjboxStats(x, coef = coef, a = a, b = b,
                          do.conf = FALSE, do.out = FALSE)$fence[2]
}

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
#' If you don't have any information about SNPs, you can try using
#'   - `infos.chr = rep(1, ncol(G))`,
#'   - `size = ncol(G)` (if SNPs are not sorted),
#'   - `roll.size = 0` (if SNPs are not sorted).
#'
#' Improvements will come in the future, as for example, warm starts in order
#' to make the SVD computations faster.
#'
#' @inheritParams bigsnpr-package
#' @inheritParams snp_clumping
#' @param fun.scaling A function that returns a named list of `mean` and `sd`
#'   for every column, to scale each of their elements such as followed:
#'   \deqn{\frac{X_{i,j} - mean_j}{sd_j}.} Default is `snp_scaleBinom()`.
#' @param k Number of singular vectors/values to compute. Default is `10`.
#'   **This algorithm should be used to compute a few singular vectors/values.**
#' @param roll.size Radius of rolling windows to smooth log-p-values.
#'   Default is `50`.
#' @param int.min.size Minimum number of consecutive outlier SNPs
#'   in order to be reported as long-range LD region. Default is `20`.
#' @param verbose Output some information on the iterations? Default is `TRUE`.
#' @param thr.r2 Threshold over the squared correlation between two SNPs.
#'   Default is `0.2`. Use `NA` if you want to skip the clumping step.
#'
#' @inherit bigstatsr::big_randomSVD return
#' @export
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
                        infos.pos = NULL,
                        ind.row = rows_along(G),
                        ind.col = cols_along(G),
                        fun.scaling = snp_scaleBinom(),
                        thr.r2 = 0.2,
                        size = 100 / thr.r2,
                        k = 10,
                        roll.size = 50,
                        int.min.size = 20,
                        is.size.in.bp = NULL,
                        ncores = 1,
                        verbose = TRUE) {

  check_args()
  assert_lengths(infos.chr, cols_along(G))
  if (!is.null(infos.pos)) assert_lengths(infos.pos, cols_along(G))

  if (!missing(is.size.in.bp))
    warning2("Parameter 'is.size.in.bp' is deprecated.")

  # Verbose?
  printf2 <- function(...) if (verbose) printf(...)

  # First clumping
  if (is.na(thr.r2)) {
    printf2("\nSkipping clumping.\n")
    ind.keep <- ind.col
  } else {
    printf2("\nPhase of clumping (on MAF) at r^2 > %s.. ", thr.r2)
    ind.keep <- snp_clumping(G, infos.chr,
                             ind.row = ind.row,
                             exclude = setdiff(cols_along(G), ind.col),
                             thr.r2 = thr.r2,
                             size = size,
                             infos.pos = infos.pos,
                             ncores = ncores)
    printf2("keep %d SNPs.\n", length(ind.keep))
  }

  iter <- 0L
  LRLDR <- LD.wiki34[0, 1:3]
  repeat {
    printf2("\nIteration %d:\n", iter <- iter + 1L)
    printf2("Computing SVD..\n")
    # SVD
    obj.svd <- big_randomSVD(G,
                             fun.scaling = fun.scaling,
                             ind.row = ind.row,
                             ind.col = ind.keep,
                             k = k,
                             ncores = ncores)

    # -log10(p-values) of being an outlier
    lpval <- -stats::predict(
      snp_pcadapt(G, obj.svd$u, ind.row, ind.keep, ncores = ncores))
    # Bonferroni-corrected threshold
    lim <- -log10(0.05 / length(lpval))

    # Roll mean to get only consecutive outliers
    ind.col.excl <- which(rollMean(lpval, size = roll.size) > lim)
    printf2("%d outlier%s detected..\n", length(ind.col.excl),
            `if`(length(ind.col.excl) > 1, "s", ""))

    # Stop or continue?
    if (length(ind.col.excl) > 0) {
      ind.keep <- ind.keep[-ind.col.excl]

      # Detection of long-range LD regions
      if (!is.null(infos.pos)) {
        # Regroup them by intervals to return long-range LD regions
        ind.range <- getIntervals(ind.col.excl, n = int.min.size)
        printf2("%d long-range LD region%s detected..\n", nrow(ind.range),
                `if`(nrow(ind.range) > 1, "s", ""))
        if (nrow(ind.range) > 0) {
          LRLDR.add <- cbind(
            infos.chr[ind.keep[ind.range[, 1]]],
            matrix(infos.pos[ind.keep[ind.range]], ncol = 2)
          )
          LRLDR[nrow(LRLDR) + rows_along(LRLDR.add), ] <- LRLDR.add
        }
      }
    } else {
      printf2("\nConverged!\n")
      break
    }
  }

  structure(obj.svd, subset = ind.keep, lrldr = LRLDR)
}

################################################################################

tukeyMC <- function(stat, coef = 1.5) {
  robustbase::adjboxStats(stat, coef = coef)$fence[2]
}


#' @rdname snp_autoSVD
#' @export
bed_autoSVD <- function(obj.bed,
                        ind.row = rows_along(obj.bed),
                        ind.col = cols_along(obj.bed),
                        thr.r2 = 0.2,
                        size = 100 / thr.r2,
                        k = 10,
                        roll.size = 50,
                        int.min.size = 20,
                        tukey.coef.row = 1.5,
                        tukey.coef.col = 3,
                        min.mac = 10,
                        ncores = 1,
                        verbose = TRUE) {

  infos.chr <- obj.bed$map$chromosome
  infos.pos <- obj.bed$map$physical.pos

  check_args()

  if (min.mac > 0) { ## TODO: option in svd instead
    stats <- bed_stats(obj.bed, ind.row, ind.col)
    mac.nok <- (stats$sum < min.mac)
    if (sum(mac.nok) > 0) {
      warning2("Discarding %d variants with MAC < %d..", sum(mac.nok), min.mac)
      ind.col <- ind.col[!mac.nok]
    }
  }

  # Verbose?
  printf2 <- function(...) if (verbose) printf(...)

  # First clumping
  if (is.na(thr.r2)) {
    printf2("\nSkipping clumping.\n")
    ind.keep <- ind.col
  } else {
    printf2("\nPhase of clumping (on MAF) at r^2 > %s.. ", thr.r2)
    ind.keep <- bed_clumping(obj.bed,
                             ind.row = ind.row,
                             exclude = setdiff(cols_along(obj.bed), ind.col),
                             thr.r2 = thr.r2,
                             size = size,
                             ncores = ncores)
    printf2("keep %d SNPs.\n", length(ind.keep))
  }

  iter <- 0L
  LRLDR <- LD.wiki34[0, 1:3]
  repeat {
    printf2("\nIteration %d:\n", iter <- iter + 1L)
    printf2("Computing SVD..\n")
    # SVD
    obj.svd <- bed_randomSVD(obj.bed,
                             ind.row = ind.row,
                             ind.col = ind.keep,
                             k = k,
                             ncores = ncores)

    # check for outlier samples
    stat <- log(dbscan::lof(obj.svd$u, k = 10))
    ind.row.excl <- which(stat > tukeyMC(stat, tukey.coef.row))
    printf2("%d outlier sample%s detected..\n", length(ind.row.excl),
            `if`(length(ind.row.excl) > 1, "s", ""))

    # -log10(p-values) of being an outlier variant
    lpval <- -stats::predict(
      bed_pcadapt(obj.bed, obj.svd$u, ind.row, ind.keep, ncores = ncores))

    # Roll mean to get only consecutive outliers
    lpval2 <- rollMean(lpval, size = roll.size)
    ind.col.excl <- which(lpval2 > tukeyMC(lpval2, tukey.coef.col))
    printf2("%d outlier variant%s detected..\n", length(ind.col.excl),
            `if`(length(ind.col.excl) > 1, "s", ""))

    # Stop or continue?
    cont <- FALSE
    if (length(ind.row.excl) > 0) {
      ind.row <- ind.row[-ind.row.excl]
      cont <- TRUE
    }
    if (length(ind.col.excl) > 0) {
      ind.keep <- ind.keep[-ind.col.excl]
      cont <- TRUE

      # Detection of long-range LD regions
      if (!is.null(infos.pos)) {
        # Regroup them by intervals to return long-range LD regions
        ind.range <- getIntervals(ind.col.excl, n = int.min.size)
        printf2("%d long-range LD region%s detected..\n", nrow(ind.range),
                `if`(nrow(ind.range) > 1, "s", ""))
        if (nrow(ind.range) > 0) {
          LRLDR.add <- cbind(
            infos.chr[ind.keep[ind.range[, 1]]],
            matrix(infos.pos[ind.keep[ind.range]], ncol = 2)
          )
          LRLDR[nrow(LRLDR) + rows_along(LRLDR.add), ] <- LRLDR.add
        }
      }
    }

    if (!cont) {
      printf2("\nConverged!\n")
      break
    }
  }

  structure(obj.svd, subset.row = ind.row, subset = ind.keep, lrldr = LRLDR)
}

################################################################################
