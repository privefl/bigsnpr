################################################################################

#' Compute quick association statistics from BGEN files
#'
#' **THIS FUNCTION WILL BE MODIFIED SOON.**
#'
#' @inheritParams snp_readBGEN
#' @param ind_row A vector of the row indices (individuals) that are used.
#'   Missing values in either `ind_row` or `y_row` are removed.
#'   **Make sure to use indices corresponding to your training set only.**
#' @param y_row A vector corresponding to `ind_row` and representing the trait
#'   with which to compute correlations.
#'   Missing values in either `ind_row` or `y_row` are removed.
#'
#' @return A list of vectors of log10(p-values) corresponding to the statistic
#'   \eqn{n \times r^2}, where r is the correlation of each variant with `y_row`.
#'
#' @import foreach
#'
#' @export
snp_assocBGEN <- function(bgenfiles, list_snp_id, y_row, ind_row,
                          bgi_dir = dirname(bgenfiles),
                          ncores = 1) {

  # Check extension of files
  sapply(bgenfiles, assert_ext, ext = "bgen")
  # Check if all files exist
  bgifiles <- file.path(bgi_dir, paste0(basename(bgenfiles), ".bgi"))
  sapply(c(bgenfiles, bgifiles), assert_exist)

  # Check list_snp_id
  assert_class(list_snp_id, "list")
  sapply(list_snp_id, assert_nona)
  assert_lengths(list_snp_id, bgenfiles)

  # Samples
  assert_lengths(y_row, ind_row)
  na_row <- is.na(ind_row) | is.na(y_row)
  if (any(na_row)) {
    printf("%d individuals removed due to missing values (out of %d).",
           sum(na_row), length(na_row))
    y_row <- y_row[!na_row]
    ind_row <- ind_row[!na_row]
  }
  if (any(duplicated(ind_row))) stop2("Can't have duplicated samples.")
  y_row <- y_row[order(ind_row)]
  N <- readBin(bgenfiles[1], what = 1L, size = 4, n = 4)[4]

  if (ncores == 1) {
    registerDoSEQ()
  } else {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }
  r2 <- foreach(ic = seq_along(bgenfiles)) %dopar% {

    snp_id <- format_snp_id(list_snp_id[[ic]])
    infos <- snp_readBGI(bgifiles[ic], snp_id)

    # Get r2
    r2 <- r2_bgen(
      filename = bgenfiles[ic],
      offsets = as.double(infos$file_start_position),
      use_ind = seq_len(N) %in% ind_row,
      decode = 0:510 / 255,
      y = y_row
    )
    r2[match(snp_id, infos$myid)]
  }

  lapply(r2, function(r2) {
    deno_y <- sum(y_row^2) - sum(y_row)^2 / length(y_row)
    stats::pchisq(length(y_row) * r2 / deno_y, df = 1, lower.tail = FALSE,
                  log.p = TRUE) / log(10)
  })
}

################################################################################
