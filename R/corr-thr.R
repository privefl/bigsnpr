################################################################################

snp_corInd <- function(Gna, list_ind,
                       ind.row = rows_along(Gna),
                       ind.col = cols_along(Gna),
                       ncores = 1) {

  assert_lengths(list_ind, ind.col)

  list_ind <- lapply(seq_along(list_ind), function(k) {
    ind <- list_ind[[k]]
    unique(sort(ind[k >= ind])) - 1L
  })

  corr <- new("dsCMatrix", uplo = "U")
  m <- length(list_ind)
  corr@Dim <- c(m, m)
  corr@i <- unlist(list_ind)
  corr@p <- c(0L, cumsum(lengths(list_ind)))
  rm(list_ind); gc()

  corr@x <- corMatInd(
    obj    = Gna,
    rowInd = ind.row,
    colInd = ind.col,
    P      = corr@p,
    I      = corr@i,
    ncores = ncores
  )

  if (anyNA(corr@x))
    warning2("NA or NaN values in the resulting correlation matrix.")

  corr
}

################################################################################

#' Extended thresholded correlation matrix
#'
#' Get all pairwise correlations `r(i,j)` for which there exists k such that
#' `r2(i,k) >= thr_r2` and `r2(j,k) >= thr_r2`. For other indices, 0s are returned.
#'
#' @inheritParams bigsnpr-package
#' @param size For one SNP, window size around this SNP to compute correlations.
#' This is internally multiplied by 1000 (e.g. `500` means 500 kb).
#' @param thr_r2 Threshold to apply on squared correlations.
#'
#' @return A sparse symmetric matrix.
#'
#' @import Matrix
#'
#' @examples
#' test <- snp_attachExtdata()
#' POS <- test$map$physical.pos
#' G <- test$genotypes
#' cor(G[, 1:5])
#' snp_cor_extendedThr(G, thr_r2 = 0.002, ind.col = 1:5,
#'                   infos.pos = POS[1:5], size = 500)
#'
#' @export
#'
snp_cor_extendedThr <- function(Gna, thr_r2, infos.pos, size,
                                ind.row = rows_along(Gna),
                                ind.col = cols_along(Gna),
                                ncores = 1) {

  message2("Initial computation of all pairwise r2 > %s..", thr_r2)
  corr <- snp_cor(Gna, thr_r2 = thr_r2, infos.pos = infos.pos, size = size,
                  ind.row = ind.row, ind.col = ind.col, ncores = ncores)

  ind <- Matrix::which(corr != 0, arr.ind = TRUE)
  prev_list_keep <- split(ind[, "row"], factor(cols_along(corr))[ind[, "col"]])
  rm(ind)

  repeat {
    corr_T <- corr %>%
      Matrix::drop0(tol = sqrt(thr_r2)) %>%
      methods::as("generalMatrix")
    new_list_keep <- find_indirect_corr(corr_T@p, corr_T@i, ncores = ncores)

    for (k in seq_along(prev_list_keep)) {
      diff_ind <- setdiff(new_list_keep[[k]], prev_list_keep[[k]])
      prev_list_keep[[k]] <- new_list_keep[[k]]
      new_list_keep[[k]] <- diff_ind
    }

    nb_to_add <- sum(lengths(new_list_keep))
    if (nb_to_add == 0) break

    message2("Computing %s additional pairwise correlations..", nb_to_add)
    corr <- corr + snp_corInd(Gna, list_ind = new_list_keep, ncores = ncores,
                              ind.row = ind.row, ind.col = ind.col)
  }

  corr
}

################################################################################

################################################################################
