################################################################################

compute_cost <- function(block_num, corr.tril, thr_r2) {
  corr.tril %>%
    as("dgTMatrix") %>%
    { .@x[block_num[.@i + 1L] != block_num[.@j + 1L]]^2 } %>%
    { sum(.[. >= thr_r2]) }
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
#' @param grid_param Grid of parameters to consider, a data frame with columns
#'   - `$min_size`: Minimum number of variants in each block.
#'   - `$max_size`: Maximum number of variants in each block.
#'   - `$lambda`: Penalty coefficient to apply on the size of the blocks.
#'     Using `0` would disable this. You can try multiple values, e.g.
#'     `c(0, 0.001, 0.01, 0.1)`.
#'
#' @return Input `grid_param` as an ordered tibble with six extra columns:
#'   - `$n_block`: Number of blocks.
#'   - `$cost`: The sum of squared correlations outside the blocks.
#'   - `$block_num`: Resulting block numbers for each variant.
#'   - `$all_last`: Last index of each block.
#'   - `$all_size`: Sizes of the blocks.
#'   - `$all_cost`: Internal costs computed.
#' This is ordered by minimum cost and maximum number of blocks.
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @example examples/example-split-LD.R
#'
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
      # L now has an extra column with all 0s for convenience
      { Matrix::sparseMatrix(i = .$i, j = .$j, x = .$x, dims = c(m, m + 1),
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

    all_last <- unlist(all_last)
    block_num <- rowSums(outer(cols_along(corr), all_last, ">")) + 1L
    all_size <- diff(c(0, all_last))
    stopifnot(all(all_size >= MIN_M & all_size <= MAX_M))

    grid_param$n_block[ic]   <- length(all_last)
    grid_param$cost[ic]      <- compute_cost(block_num, corr, thr_r2)
    grid_param$block_num[ic] <- list(block_num)
    grid_param$all_last[ic]  <- list(all_last)
    grid_param$all_size[ic]  <- list(all_size)
    grid_param$all_cost[ic]  <- list(cost_path$C)
  }

  tibble::as_tibble(grid_param[order(grid_param$cost, -grid_param$n_block), ])
}

################################################################################
