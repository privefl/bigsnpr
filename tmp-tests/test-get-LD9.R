corr <- readRDS(url("https://www.dropbox.com/s/65u96jf7y32j2mj/spMat.rds?raw=1"))
m <- ncol(corr)
POS <- seq_len(m)
min_size <- 9
max_size <- 49
max_K <- 50
thr_r2 <- 0.01

corr <- Matrix::tril(corr)

Rcpp::sourceCpp('src/split-LD.cpp')

# Precomputing L, E and computing all cost paths
library(magrittr)
cost_path <- corr %>%
  { get_L(.@p, .@i, .@x, thr_r2 = thr_r2) } %>%
  # L now has an extra column with all 0s for convenience
  { Matrix::sparseMatrix(i = .$i, j = .$j, x = .$x, dims = c(m, m + 1),
                         triangular = FALSE, index1 = FALSE) } %>%
  get_C(pos = POS, min_size = min_size, max_size = max_size, K = max_K)

# Reconstructing paths
res <- do.call("rbind", lapply(1:max_K, function(K) {

  print(K)

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

  all_first <- c(1L, head(all_last + 1L, -1))
  all_size <- POS[all_last] - POS[all_first]
  stopifnot(all(all_size >= min_size & all_size <= max_size))

  block_num <- rowSums(outer(1:m, all_last, ">")) + 1L

  tibble::tibble(
    n_block   = length(all_last),
    cost      = cost,
    block_num = list(block_num),
    all_last  = list(all_last),
    all_size  = list(all_size)
  )
}))


res2 <- dplyr::bind_cols(k = 9:40, res)
ggplot2::qplot(k, n_block, data = res2)


ggplot2::qplot(n_block, cost, data = res) + ggplot2::scale_y_log10()
