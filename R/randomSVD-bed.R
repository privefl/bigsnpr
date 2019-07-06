################################################################################

# Parallel implementation
svds4.par <- function(obj.bed, ind.row, ind.col, k, tol, verbose, ncores) {

  n <- length(ind.row)
  m <- length(ind.col)
  intervals <- CutBySize(m, nb = ncores)

  TIME <- 0.001

  Ax   <- FBM(n, ncores)
  Atx  <- FBM(m, 1)
  calc <- FBM(ncores, 1, init = 0)
  nb_nona_row <- FBM(n, ncores)

  cluster_type <- getOption("bigstatsr.cluster.type")
  if (verbose) {
    cl <- parallel::makeCluster(1 + ncores, type = cluster_type, outfile = "")
  } else {
    cl <- parallel::makeCluster(1 + ncores, type = cluster_type)
  }
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  res <- foreach(ic = 0:ncores) %dopar% {

    if (ic == 0) { # I'm the master

      printf <- function(...) cat(sprintf(...))
      it <- 0
      # A
      A <- function(x, args) {
        printf("%d - computing A * x\n", it <<- it + 1)
        Atx[] <- x
        calc[] <- 1  # make them work
        # Master wait for its slaves to finish working
        while (sum(calc[]) > 0) Sys.sleep(TIME)
        rowSums(Ax[]) * m / rowSums(nb_nona_row[])
      }
      # Atrans
      Atrans <- function(x, args) {
        printf("%d - computing At * x\n", it <<- it + 1)
        Ax[, 1] <- x
        calc[] <- 2  # make them work
        # Master wait for its slaves to finish working
        while (sum(calc[]) > 0) Sys.sleep(TIME)
        Atx[]
      }

      res <- RSpectra::svds(A, k, nu = k, nv = k, opts = list(tol = tol),
                            Atrans = Atrans, dim = c(n, m))

      calc[] <- 3 # end

      res
    } else { # You're my slaves
      # Get their part
      lo <- intervals[ic, "lower"]
      up <- intervals[ic, "upper"]
      ind.col.part <- ind.col[lo:up]

      # Scaling
      stats <- bigsnpr:::bed_stats(obj.bed, ind.row, ind.col.part)
      af <- stats$sum / (2 * stats$nb_nona_col)
      center <- 2 * af
      scale <- sqrt(2 * af * (1 - af))
      nb_nona_row[, ic] <- stats$nb_nona_row

      repeat {
        # Slaves wait for their master to give them orders
        while (calc[ic] == 0) Sys.sleep(TIME)
        c <- calc[ic]
        # Slaves do the hard work
        if (c == 1) {
          # Compute A * x (part)
          Ax[, ic] <- pMatVec4(obj.bed, ind.row, ind.col.part, center, scale,
                               Atx[lo:up]) #* m / stats$nb_nona_row
        } else if (c == 2) {
          # Compute At * x (part)
          Atx[lo:up] <- cpMatVec4(obj.bed, ind.row, ind.col.part, center, scale,
                                  Ax[, 1]) * n / stats$nb_nona_col
        } else if (c == 3) {
          # End
          break
        } else {
          stop("RandomSVD: unclear order from the master.")
        }
        calc[ic] <- 0
      }

      data.frame(center = center, scale = scale)
    }
  }

  # Separate the results and combine the scaling vectors
  l <- do.call("c", res[-1])
  res <- res[[1]]
  s <- c(TRUE, FALSE)
  res$center <- unlist(l[s],  use.names = FALSE)
  res$scale  <- unlist(l[!s], use.names = FALSE)

  # Return
  res
}

################################################################################

# Single core implementation
svds4.seq <- function(obj.bed, ind.row, ind.col, k, tol, verbose) {

  n <- length(ind.row)
  m <- length(ind.col)

  # scaling
  stats <- bigsnpr:::bed_stats(obj.bed, ind.row, ind.col)
  af <- stats$sum / (2 * stats$nb_nona_col)
  center <- 2 * af
  scale <- sqrt(2 * af * (1 - af))

  printf <- function(...) if (verbose) cat(sprintf(...))
  it <- 0
  # A
  A <- function(x, args) {
    printf("%d - computing A * x\n", it <<- it + 1)
    pMatVec4(obj.bed, ind.row, ind.col, center, scale, x) *
      m / stats$nb_nona_row
  }
  # Atrans
  Atrans <- function(x, args) {
    printf("%d - computing At * x\n", it <<- it + 1)
    cpMatVec4(obj.bed, ind.row, ind.col, center, scale, x) *
      n / stats$nb_nona_col
  }

  res <- RSpectra::svds(A, k, nu = k, nv = k, opts = list(tol = tol),
                        Atrans = Atrans, dim = c(n, m))

  res$center <- center
  res$scale  <- scale

  res
}

################################################################################

#' Randomized partial SVD
#'
#' An algorithm for partial SVD (or PCA) of a Filebacked Big Matrix based on the
#' algorithm in RSpectra (by Yixuan Qiu and Jiali Mei).\cr
#' This algorithm is linear in time in all dimensions and is very
#' memory-efficient. Thus, it can be used on very large bed files.
#'
#' @inheritParams bigstatsr-package
#' @param k Number of singular vectors/values to compute. Default is `10`.
#' **This algorithm should be used to compute only a
#' few singular vectors/values.**
#' @param tol Precision parameter of [svds][RSpectra::svds].
#' Default is `1e-4`.
#' @param verbose Should some progress be printed? Default is `FALSE`.
#'
#' @export
#'
#' @return A named list (an S3 class "big_SVD") of
#' - `d`, the singular values,
#' - `u`, the left singular vectors,
#' - `v`, the right singular vectors,
#' - `niter`, the number of the iteration of the algorithm,
#' - `nops`, number of Matrix-Vector multiplications used,
#' - `center`, the centering vector,
#' - `scale`, the scaling vector.
#'
#' Note that to obtain the Principal Components, you must use
#' [predict][predict.big_SVD] on the result. See examples.
#'
#' @example examples/example-randomSVD.R
#' @seealso [svds][RSpectra::svds]
#'
bed_randomSVD <- function(
  obj.bed,
  ind.row = rows_along(obj.bed),
  ind.col = cols_along(obj.bed),
  k = 10,
  tol = 1e-4,
  verbose = FALSE,
  ncores = 1
) {

  check_args()

  if (ncores > 1) {
    res <- svds4.par(obj.bed, ind.row, ind.col, k, tol, verbose, ncores)
  } else {
    res <- svds4.seq(obj.bed, ind.row, ind.col, k, tol, verbose)
  }

  structure(res, class = "big_SVD")
}

################################################################################
