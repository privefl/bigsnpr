################################################################################

verif_cost <- function(cost_to_verif, block_num, corr.tril, thr_r2) {

  cost <- corr.tril %>%
    as("dgTMatrix") %>%
    { .@x[block_num[.@i + 1L] != block_num[.@j + 1L]]^2 } %>%
    { sum(.[. >= thr_r2]) }

  if (!isTRUE(all.equal(cost, cost_to_verif, tolerance = 1e-6)))
    stop2("Computed costs are not matching (%s != %s).\n  %s",
          cost, cost_to_verif, "Please report this issue.")
}

################################################################################

#' Independent LD blocks
#'
#' Split a correlation matrix in blocks as independent as possible.
#' This will find the splitting in blocks that minimize the sum of squared
#' correlation between these blocks (i.e. everything outside these blocks).
#'
#' @param corr Sparse correlation matrix. Usually the output of [snp_cor()].
#' @param thr_r2 Threshold under which squared correlations are ignored.
#' @param grid_param Grid of parameters to consider, a data frame with
#'   - `$min_size`: Minimum number of variants in each block.
#'   - `$max_size`: Maximum number of variants in each block.
#'   - `$lambda`: Penalty coefficient to apply on the size of the blocks.
#'     Using `0` would disable this, while using `1e-3` could be a good value
#'     to get more smaller blocks.
#'
#' @importFrom magrittr %>%
#'
#' @return Input `grid_param` as tibble with five extra columns:
#'   - `$all_cost`: Internal costs computed.
#'   - `$cost`: The first cost, `C[1]`, represents the sum of squared
#'     correlations outside the blocks.
#'   - `$all_last`: Last index of each block.
#'   - `$all_size`: Sizes of the blocks.
#'   - `$n_block`: Number of blocks.
#' @export
#'
#' @examples
snp_ldsplit <- function(corr, thr_r2, grid_param) {

  assert_df_with_names(grid_param, c("min_size", "max_size", "lambda"))

  m <- ncol(corr)
  corr <- Matrix::tril(corr)

  for (ic in rows_along(grid_param)) {

    MIN_M  <- grid_param$min_size[ic]
    MAX_M  <- grid_param$max_size[ic]
    LAMBDA <- grid_param$lambda[ic]

    # Precomputing L, E and computing all cost paths
    cost_path <- corr %>%
      { get_L(.@p, .@i, .@x, thr_r2 = thr_r2) } %>%
      { Matrix::sparseMatrix(i = .$i, j = .$j, x = .$x, dims = c(m, m),
                             triangular = FALSE, index1 = FALSE) } %>%
      get_C(min_size = MIN_M, max_size = MAX_M, lambda = LAMBDA)

    # Reconstruct path
    all_last <- list()
    j <- 0
    repeat {
      j <- cost_path$best_ind[j + 1L]
      if (is.na(j)) break
      all_last[[length(all_last) + 1L]] <- j
    }
    all_last <- c(unlist(all_last), m)

    block_num <- rowSums(outer(cols_along(corr), all_last, ">")) + 1L
    verif_cost(cost_path$C[1], block_num, corr, thr_r2)

    all_size <- diff(c(0, all_last))
    stopifnot(all(all_size >= MIN_M & all_size <= MAX_M))

    grid_param$n_block[ic]   <- length(all_last)
    grid_param$cost[ic]      <- cost_path$C[1]
    grid_param$block_num[ic] <- list(block_num)
    grid_param$all_last[ic]  <- list(all_last)
    grid_param$all_size[ic]  <- list(all_size)
    grid_param$all_cost[ic]  <- list(cost_path$C)
  }

  tibble::as_tibble(grid_param)
}

################################################################################
