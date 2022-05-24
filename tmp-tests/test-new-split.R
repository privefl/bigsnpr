library(bigsnpr)
bedfile <- download_1000G("../datasets")
obj.bed <- bed(bedfile)
fam2 <- bigreadr::fread2(sub_bed(bedfile, ".fam2"))
ind <- which(fam2$`Super Population` == "EUR")
ind_chr <- which(obj.bed$map$chromosome == 4)

corr0 <- bed_cor(obj.bed, ind.row = ind, ind.col = ind_chr, ncores = nb_cores())
corr <- Matrix::tril(corr0)
# corr@x <- ifelse(is.na(corr@x), 0, corr@x)
# corr@x <- ifelse(corr@x^2 < 0.05, 0, corr@x)

L <- bigsnpr:::get_L(corr@p, corr@i, corr@x, thr_r2 = 0.05)
anyNA(L$x)
length(L$x) / length(corr@x)

min_size <- 200
max_size <- 10000
max_K <- 80

reconstruct_paths <- function(cost_path) {

  do.call("rbind", lapply(1:max_K, function(K) {

    cost <- cost_path$C[1, K]
    if (is.na(cost) || is.infinite(cost)) return(NULL)

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
      perc_kept = bigsnpr:::get_perc(corr@p, corr@i, block_num),
      block_num = list(block_num),
      all_last  = list(all_last),
      all_size  = list(all_size)
    )
  }))
}

# X <- FBM(1, 1, "float")
# X[1] <- Inf
# X[]

library(magrittr)
m <- ncol(corr)
system.time(
  cost_path <- corr %>%
    { bigsnpr:::get_L(.@p, .@i, .@x, thr_r2 = 0.05) } %>%
    # L now has an extra column with all 0s for convenience
    { Matrix::sparseMatrix(i = .$i, j = .$j, x = .$x, dims = c(m, m + 1),
                           triangular = FALSE, index1 = FALSE) } %>%
    bigsnpr:::get_C(min_size = min_size, max_size = max_size, K = max_K)
)
# 10 / 64 / 148 sec

Rcpp::sourceCpp('tmp-tests/test-new-split2.cpp')

system.time(
  cost_path2 <- corr %>%
    { get_L2(.@p, .@i, .@x, thr_r2 = 0.05, max_r2 = 0.5) } %>%
    # L now has an extra column with all 0s for convenience
    { Matrix::sparseMatrix(i = .$i, j = .$j, x = .$x, dims = c(m, m + 1),
                           triangular = FALSE, index1 = FALSE) } %>%
    get_C2(min_size = min_size, max_size = max_size, max_K = max_K)
)
# 0.7 / 2 / 3 sec


mean(is.na(cost_path$C))
mean(is.infinite(cost_path2$C))
# plot(cost_path$C, cost_path2$C)
cost_path$C[1:5, 1:5]
cost_path2$C[1:5, 1:5]



library(ggplot2)
library(dplyr)

all_splits <- bind_rows(
  bind_cols(Method = "Normal", reconstruct_paths(cost_path)),
  bind_cols(Method = "With max r2", reconstruct_paths(cost_path2))
)

qplot(data = all_splits, n_block, cost, color = Method, size = Method) +
  theme_bw(12) +
  theme(legend.position = "top") + guides(colour = guide_legend(nrow = 1)) +
  scale_size_manual(values = c(1, 0.5, 0.5) * 2) +
  ylim(0, 500) +
  labs(y = "Sum of squared correlations outside blocks (log-scale)",
       x = "Number of blocks")

qplot(data = all_splits, perc_kept, cost, color = Method, size = Method) +
  theme_bw(12) +
  theme(legend.position = "top") + guides(colour = guide_legend(nrow = 1)) +
  scale_size_manual(values = c(1, 0.5, 0.5) * 2) +
  labs(y = "Sum of squared correlations outside blocks (log-scale)",
       x = "% of non-zero values kept") +
  # ylim(0, min(quantile(all_splits$cost, 0.7), 500)) +
  geom_vline(xintercept = 0.6, linetype = 3) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2))


compute_cost <- function(block_num, corr.tril, thr_r2) {
  corr.tril %>%
    methods::as("dgTMatrix") %>%
    { .@x[block_num[.@i + 1L] != block_num[.@j + 1L]]^2 } %>%
    { sum(.[. >= thr_r2]) }
}
(res <- filter(all_splits, n_block == 30))
compute_cost(res$block_num[[1]], corr, 0.05)
compute_cost(res$block_num[[2]], corr, 0.05)
