################################################################################

compute_cost <- function(block_num, corr.tril, thr_r2) {
  corr.tril %>%
    methods::as("dgTMatrix") %>%
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
#' @param corr Sparse correlation matrix. Usually, the output of [snp_cor()].
#' @param thr_r2 Threshold under which squared correlations are ignored.
#'   This is useful to avoid counting noise, which should give clearer patterns
#'   of costs vs. number of blocks. It is therefore possible to have a splitting
#'   cost of 0. If this parameter is used, then `corr` can be computed using the
#'   same parameter in [snp_cor()] (to increase the sparsity of the resulting matrix).
#' @param min_size Minimum number of variants in each block. This is used not to
#'   have a disproportionate number of small blocks.
#' @param max_size Maximum number of variants in each block. This is used not to
#'   have blocks that are too large, e.g. to limit computational and memory
#'   requirements of applications that would use these blocks. For some long-range
#'   LD regions, it may be needed to allow for large blocks.
#' @param max_K Maximum number of blocks to consider. All optimal solutions for K
#'   from 1 to `max_K` will be returned. Some of these K might not have any corresponding
#'   solution due to the limitations in size of the blocks. For example, splitting
#'   10,000 variants in blocks with at least 500 and at most 2000 variants implies
#'   that there are at least 5 and at most 20 blocks. Then, the choice of K depends
#'   on the application, but a simple solution is to choose the largest K for which
#'   the cost is lower than some threshold.
#'
#' @return A tibble with five columns:
#'   - `$n_block`: Number of blocks.
#'   - `$cost`: The sum of squared correlations outside the blocks.
#'   - `$perc_kept`: Percentage of initial non-zero values kept within the blocks defined.
#'   - `$block_num`: Resulting block numbers for each variant.
#'   - `$all_last`: Last index of each block.
#'   - `$all_size`: Sizes of the blocks.
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @example examples/example-split-LD.R
#'
snp_ldsplit <- function(corr, thr_r2, min_size, max_size, max_K) {

  m <- ncol(corr)
  stopifnot(min_size >= 1 && max_size <= m)

  corr <- Matrix::tril(corr)

  # Precomputing L, E and computing all cost paths
  cost_path <- corr %>%
    { get_L(.@p, .@i, .@x, thr_r2 = thr_r2) } %>%
    # L now has an extra column with all 0s for convenience
    { Matrix::sparseMatrix(i = .$i, j = .$j, x = .$x, dims = c(m, m + 1),
                           triangular = FALSE, index1 = FALSE) } %>%
    get_C(min_size = min_size, max_size = max_size, K = max_K)

  # Reconstructing paths
  do.call("rbind", lapply(1:max_K, function(K) {

    cost <- cost_path$C[1, K]
    if (is.na(cost)) return(NULL)

    all_last <- list()
    j <- 0
    k <- K
    repeat {
      j <- cost_path$best_ind[j + 1L, k]
      all_last[[length(all_last) + 1L]] <- j
      if (k == 1) break
      k <- k - 1L
    }

    all_last <- unlist(all_last)
    stopifnot(length(all_last) == K)

    all_size <- diff(c(0, all_last))
    stopifnot(all(all_size >= min_size & all_size <= max_size))

    block_num <- rowSums(outer(1:m, all_last, ">")) + 1L

    tibble::tibble(
      n_block   = length(all_last),
      cost      = cost,
      perc_kept = get_perc(corr@p, corr@i, block_num),
      block_num = list(block_num),
      all_last  = list(all_last),
      all_size  = list(all_size)
    )
  }))
}

################################################################################
