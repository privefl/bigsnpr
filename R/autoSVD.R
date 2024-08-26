################################################################################

# regroup consecutive integers in intervals
getIntervals <- function(x, n = 2) {

  le <- rle(diff(x))

  ind <- which((le$values == 1) & (le$lengths >= (n - 1)))
  pos <- cumsum(le$lengths) + 1

  cbind(x[c(1, pos)[ind]], x[pos[ind]])
}

################################################################################

#' Truncated SVD while limiting LD
#'
#' Fast truncated SVD with initial pruning and that iteratively removes
#' long-range LD regions. Some variants are removing due to the initial clumping,
#' then more and more variants are removed at each iteration. You can access the
#' indices of the remaining variants with `attr(*, "subset")`. If some of the
#' variants removed are contiguous, the regions are reported in `attr(*, "lrldr")`.
#'
#' If you don't have any information about variants, you can try using
#'   - `infos.chr = rep(1, ncol(G))`,
#'   - `size = ncol(G)` (if variants are not sorted),
#'   - `roll.size = 0` (if variants are not sorted).
#'
#' @inheritParams bigsnpr-package
#' @inheritParams snp_clumping
#' @param fun.scaling A function with parameters `X` (or `obj.bed`), `ind.row` and
#'   `ind.col`, and that returns a data.frame with `$center` and `$scale` for the
#'   columns corresponding to `ind.col`, to scale each of their elements such as followed:
#'   \deqn{\frac{X_{i,j} - center_j}{scale_j}.} Default uses binomial scaling.
#'   You can also provide your own `center` and `scale` by using [as_scaling_fun()].
#' @param k Number of singular vectors/values to compute. Default is `10`.
#'   **This algorithm should be used to compute a few singular vectors/values.**
#' @param roll.size Radius of rolling windows to smooth log-p-values.
#'   Default is `50`.
#' @param int.min.size Minimum number of consecutive outlier variants
#'   in order to be reported as long-range LD region. Default is `20`.
#' @param verbose Output some information on the iterations? Default is `TRUE`.
#' @param thr.r2 Threshold over the squared correlation between two variants.
#'   Default is `0.2`. Use `NA` if you want to skip the clumping step.
#' @param alpha.tukey Default is `0.1`. The type-I error rate in outlier
#'   detection (that is further corrected for multiple testing).
#' @param min.mac Minimum minor allele count (MAC) for variants to be included.
#'   Default is `10`.
#' @param max.iter Maximum number of iterations of outlier detection.
#'   Default is `5`.
#'
#' @inherit bigstatsr::big_randomSVD return
#' @export
#'
#' @examples
#' ex <- snp_attachExtdata()
#' G <- ex$genotypes
#'
#' obj.svd <- snp_autoSVD(G,
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
                        alpha.tukey = 0.05,
                        min.mac = 10,
                        max.iter = 5,
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

  if (min.mac > 0) {
    maf <- snp_MAF(G, ind.row, ind.col, ncores = ncores)
    min.maf <- min.mac / (2 * length(ind.row))
    mac.nok <- (maf < min.maf)
    printf2("Discarding %d variant%s with MAC < %s.\n", sum(mac.nok),
            `if`(sum(mac.nok) > 1, "s", ""), min.mac)
    ind.keep <- ind.col[!mac.nok]
  } else {
    stop2("You cannot use variants with no variation; set min.mac > 0.")
  }

  # First clumping
  if (is.na(thr.r2)) {
    printf2("\nSkipping clumping.\n")
  } else {
    printf2("\nPhase of clumping (on MAF) at r^2 > %s.. ", thr.r2)
    ind.keep <- snp_clumping(G, infos.chr,
                             ind.row = ind.row,
                             exclude = setdiff(cols_along(G), ind.keep),
                             thr.r2 = thr.r2,
                             size = size,
                             infos.pos = infos.pos,
                             ncores = ncores)
    printf2("keep %d variants.\n", length(ind.keep))
  }

  iter <- 0L
  LRLDR <- data.frame(Chr = integer(), Start = integer(), Stop = integer(), Iter = integer())
  repeat {
    iter <- iter + 1L
    printf2("\nIteration %d:\n", iter)
    printf2("Computing SVD..\n")
    # SVD
    obj.svd <- big_randomSVD(G,
                             fun.scaling = fun.scaling,
                             ind.row = ind.row,
                             ind.col = ind.keep,
                             k = k,
                             ncores = ncores)

    if (iter > max.iter) {
      printf2("Maximum number of iterations reached.\n")
      break
    }

    # check for outlier variants
    S.col <- sqrt(bigutilsr::dist_ogk(obj.svd$v))
    # roll mean to get only consecutive outliers (by chromosome)
    ind.split <- split(seq_along(S.col), infos.chr[ind.keep])
    S2.col <- double(length(S.col))
    for (ind in ind.split)
      S2.col[ind] <- bigutilsr::rollmean(S.col[ind], roll.size)
    S2.col.thr <- bigutilsr::tukey_mc_up(S2.col, alpha = alpha.tukey)
    ind.col.excl <- which(S2.col > S2.col.thr)
    printf2("%d outlier variant%s detected..\n", length(ind.col.excl),
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
        for (i in rows_along(ind.range)) {
          seq.range <- seq2(ind.range[i, ])
          seq.range.chr <- infos.chr[ind.keep[seq.range]]
          chr <- names(sort(table(seq.range.chr), decreasing = TRUE)[1])  ## mode
          in.chr <- (seq.range.chr == chr)
          range.in.chr <- range(infos.pos[ind.keep[seq.range[in.chr]]])
          LRLDR[nrow(LRLDR) + 1L, ] <-
            data.frame(seq.range.chr[in.chr][1], range.in.chr[1], range.in.chr[2], iter)
        }
      }
    } else {
      printf2("\nConverged!\n")
      break
    }
  }

  structure(
    obj.svd,
    subset = ind.keep,
    lrldr = LRLDR[with(LRLDR, order(Chr, Start, Stop)), ]
  )
}

################################################################################

#' Randomized partial SVD
#'
#' Partial SVD (or PCA) of a genotype matrix stored as a PLINK (.bed) file.#'
#'
#' @inheritParams bigsnpr-package
#' @inherit bigstatsr::big_randomSVD params return
#'
#' @export
#'
#' @examples
#' bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
#' obj.bed <- bed(bedfile)
#'
#' str(bed_randomSVD(obj.bed))
#'
bed_randomSVD <- function(
  obj.bed,
  fun.scaling = bed_scaleBinom,
  ind.row = rows_along(obj.bed),
  ind.col = cols_along(obj.bed),
  k = 10,
  tol = 1e-4,
  verbose = FALSE,
  ncores = 1
) {

  big_randomSVD(obj.bed$light, fun.scaling, ind.row, ind.col,
                k, tol, verbose, ncores,
                bed_prodVec, bed_cprodVec)
}

################################################################################

#' @rdname snp_autoSVD
#'
#' @export
bed_autoSVD <- function(obj.bed,
                        ind.row = rows_along(obj.bed),
                        ind.col = cols_along(obj.bed),
                        fun.scaling = bed_scaleBinom,
                        thr.r2 = 0.2,
                        size = 100 / thr.r2,
                        k = 10,
                        roll.size = 50,
                        int.min.size = 20,
                        alpha.tukey = 0.05,
                        min.mac = 10,
                        max.iter = 5,
                        ncores = 1,
                        verbose = TRUE) {

  infos.chr <- obj.bed$map$chromosome
  infos.pos <- obj.bed$map$physical.pos

  check_args()

  # Verbose?
  printf2 <- function(...) if (verbose) printf(...)

  if (min.mac > 0) {
    mac <- bed_MAF(obj.bed, ind.row, ind.col, ncores = ncores)$mac
    mac.nok <- (mac < min.mac)
    printf2("Discarding %d variant%s with MAC < %s.\n", sum(mac.nok),
            `if`(sum(mac.nok) > 1, "s", ""), min.mac)
    ind.keep <- ind.col[!mac.nok]
  } else {
    stop2("You cannot use variants with no variation; set min.mac > 0.")
  }

  # First clumping
  if (is.na(thr.r2)) {
    printf2("\nSkipping clumping.\n")
  } else {
    printf2("\nPhase of clumping (on MAC) at r^2 > %s.. ", thr.r2)
    ind.keep <- bed_clumping(obj.bed,
                             ind.row = ind.row,
                             exclude = setdiff(cols_along(obj.bed), ind.keep),
                             thr.r2 = thr.r2,
                             size = size,
                             ncores = ncores)
    printf2("keep %d variants.\n", length(ind.keep))
  }

  iter <- 0L
  LRLDR <- data.frame(Chr = integer(), Start = integer(), Stop = integer(), Iter = integer())
  repeat {
    iter <- iter + 1L
    printf2("\nIteration %d:\n", iter)
    printf2("Computing SVD..\n")
    # SVD
    obj.svd <- bed_randomSVD(obj.bed,
                             fun.scaling = fun.scaling,
                             ind.row = ind.row,
                             ind.col = ind.keep,
                             k = k,
                             ncores = ncores)

    if (iter > max.iter) {
      printf2("Maximum number of iterations reached.\n")
      break
    }

    # check for outlier variants
    S.col <- sqrt(bigutilsr::dist_ogk(obj.svd$v))
    # roll mean to get only consecutive outliers (by chromosome)
    ind.split <- split(seq_along(S.col), infos.chr[ind.keep])
    S2.col <- double(length(S.col))
    for (ind in ind.split)
      S2.col[ind] <- bigutilsr::rollmean(S.col[ind], roll.size)
    S2.col.thr <- bigutilsr::tukey_mc_up(S2.col, alpha = alpha.tukey)
    ind.col.excl <- which(S2.col > S2.col.thr)
    printf2("%d outlier variant%s detected..\n", length(ind.col.excl),
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
        for (i in rows_along(ind.range)) {
          seq.range <- seq2(ind.range[i, ])
          seq.range.chr <- infos.chr[ind.keep[seq.range]]
          chr <- names(sort(table(seq.range.chr), decreasing = TRUE)[1])  ## mode
          in.chr <- (seq.range.chr == chr)
          range.in.chr <- range(infos.pos[ind.keep[seq.range[in.chr]]])
          LRLDR[nrow(LRLDR) + 1L, ] <-
            data.frame(seq.range.chr[in.chr][1], range.in.chr[1], range.in.chr[2], iter)
        }
      }
    } else {
      printf2("\nConverged!\n")
      break
    }
  }

  structure(
    obj.svd,
    subset = ind.keep,
    lrldr = LRLDR[with(LRLDR, order(Chr, Start, Stop)), ]
  )
}

################################################################################
