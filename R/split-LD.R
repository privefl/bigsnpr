################################################################################

reconstruct_paths <- function(corr.p, corr.i, cost_path, min_size, max_size,
                              max_cost, prev_costs) {

  do.call("rbind", lapply(cols_along(cost_path$C), function(K) {

    cost <- cost_path$C[1, K]
    if (cost > max_cost) return(NULL)
    if (cost < prev_costs[K]) {
      prev_costs[K] <- cost  # change globally because use an FBM
    } else return(NULL)

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

    tibble::tibble(
      max_size  = max_size,
      n_block   = K,
      cost      = cost,
      cost2     = sum(all_size^2),
      perc_kept = get_perc(corr.p, corr.i, all_last = all_last - 1L),
      all_last  = list(all_last),
      all_size  = list(all_size)
    )
  }))
}

################################################################################

#' Independent LD blocks
#'
#' Split a correlation matrix in blocks as independent as possible.
#' This finds the splitting in blocks that minimizes the sum of squared
#' correlation between these blocks (i.e. everything outside these blocks).
#' In case of equivalent splits, it then minimizes the sum of squared sizes
#' of the blocks.
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
#'   LD regions, it may be needed to allow for large blocks. You can now provide a
#'   vector of values to try.
#' @param max_K Maximum number of blocks to consider. All optimal solutions for K
#'   from 1 to `max_K` will be returned. Some of these K might not have any corresponding
#'   solution due to the limitations in size of the blocks. For example, splitting
#'   10,000 variants in blocks with at least 500 and at most 2000 variants implies
#'   that there are at least 5 and at most 20 blocks. Then, the choice of K depends
#'   on the application, but a simple solution is to choose the largest K for which
#'   the cost is lower than some threshold. Default is `500`.
#' @param max_r2 Maximum squared correlation allowed for one pair of variants in
#'   two different blocks. This is used to make sure that strong correlations are
#'   not discarded and also to speed up the algorithm. Default is `0.3`.
#' @param max_cost Maximum cost reported. Default is `ncol(corr) / 200`.
#' @param pos_scaled Vector of positions. The positions should be scaled so that
#'   limits of a block must be separated by a distance of 1 at the maximum. E.g.
#'   if the positions are in base pairs (bp), and you want a maximum distance of
#'   10 Mbp, you need to provide the vector of positions divided by 10e6.
#'
#' @return Either `NULL` when no block splitting satisfies the conditions,
#'   or a tibble with seven columns:
#'   - `$max_size`: Input parameter, useful when providing a vector of values to try.
#'   - `$n_block`: Number of blocks.
#'   - `$cost`: The sum of squared correlations outside the blocks.
#'   - `$cost2`: The sum of squared sizes of the blocks.
#'   - `$perc_kept`: Percentage of initial non-zero values kept within the blocks defined.
#'   - `$all_last`: Last index of each block.
#'   - `$all_size`: Sizes of the blocks.
#'   - `$block_num`: Resulting block numbers for each variant. This is not reported
#'     anymore, but can be computed with `rep(seq_along(all_size), all_size)`.
#'
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @example examples/example-split-LD.R
#'
snp_ldsplit <- function(corr, thr_r2, min_size, max_size,
                        max_K = 500,
                        max_r2 = 0.3,
                        max_cost = ncol(corr) / 200,
                        pos_scaled = rep(0, ncol(corr))) {

  m <- ncol(corr)
  stopifnot(min_size >= 1 && all(max_size <= m))
  assert_lengths(pos_scaled, cols_along(corr))

  corr <- Matrix::tril(corr)
  stopifnot(all(Matrix::diag(corr) != 0))

  # to allow max_cost=Inf, and *2 to count both tri
  max_cost <- min(max_cost, crossprod(corr@x) * 2)
  prev_costs <- FBM(max_K, 1, init = Inf)

  # Precomputing L
  L <- corr %>%
    { get_L(.@p, .@i, .@x, thr_r2 = thr_r2, max_r2 = max_r2) } %>%
    # L now has an extra column with all 0s for convenience
    { Matrix::sparseMatrix(i = .$i, j = .$j, x = .$x, dims = c(m, m + 1),
                           triangular = FALSE, index1 = FALSE) }

  corr.p <- as.double(corr@p)
  corr.i <- corr@i
  rm(corr)

  do.call("rbind", lapply(sort(max_size), function(one_max_size) {

    # Precomputing E and computing all cost paths
    cost_path <- get_C(L, min_size = min_size, max_size = one_max_size,
                       max_K = max_K, max_cost = max_cost,
                       pos_scaled = pos_scaled)

    # Reconstructing paths
    reconstruct_paths(corr.p, corr.i, cost_path, min_size, one_max_size,
                      max_cost, prev_costs)
  }))
}

################################################################################
