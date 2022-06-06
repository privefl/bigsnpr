corr <- readRDS("../datasets/LD_chr4.rds")

# corr.tril <- Matrix::tril(corr)
# Rcpp::sourceCpp('src/split-LD.cpp')
#
# block_num <- sort(rep_len(1:15, ncol(corr)))
# all_last <- sapply(1:15, function(num) max(which(block_num == num)))
# get_perc(corr.tril@p, corr.tril@i, all_last = all_last - 1L)
# 0.8583042


# system.time(
#   res <- bigsnpr::snp_ldsplit(corr, thr_r2 = 0.02, min_size = 100,
#                               max_size = 8000, max_K = 80, max_r2 = 1)
# )
# 13 / 30 / 90 sec

system.time(
  res2 <- bigsnpr::snp_ldsplit(corr, thr_r2 = 0.02, min_size = 100,
                              max_size = 8000, max_K = 80, max_r2 = 0.3)
)
# 2 / 4 / 10.5 sec

system.time(
  res3 <- bigsnpr::snp_ldsplit(corr, thr_r2 = 0.02, min_size = 100,
                               max_size = 8000, max_K = 800, max_r2 = 0.3)
)

SEQ <- round(bigsnpr::seq_log(1000 + ncol(corr) / 30,
                              5000 + ncol(corr) / 10,
                              length.out = 6))
SEQ
system.time(
  res4 <- bigsnpr::snp_ldsplit(corr, thr_r2 = 0.02, min_size = 100,
                               max_size = SEQ, max_r2 = 0.1)
) # 52.6 sec -> 14 sec without get_perc
# could parallelize get_perc with OpenMP by looping over all all_last

library(ggplot2)
all_splits <- res4

qplot(data = all_splits, n_block, cost, color = as.factor(max_size)) +
  theme_bw(12) +
  scale_y_log10() +
  theme(legend.position = "top") + guides(colour = guide_legend(nrow = 1)) +
  labs(y = "Sum of squared correlations outside blocks (log-scale)",
       x = "Number of blocks", color = "Maximum block size")

qplot(data = all_splits, perc_kept, cost, color = as.factor(max_size)) +
  theme_bw(12) +
  # scale_y_log10() +
  theme(legend.position = "top") + guides(colour = guide_legend(nrow = 1)) +
  labs(y = "Sum of squared correlations outside blocks (log-scale)",
       x = "Number of blocks", color = "Maximum block size")

all_splits$cost2 <- sapply(all_splits$all_size, function(sizes) sum(sizes^2))

qplot(data = all_splits, cost2, cost, color = as.factor(max_size)) +
  theme_bw(12) +
  scale_y_log10() +
  theme(legend.position = "top") + guides(colour = guide_legend(nrow = 1)) +
  labs(y = "Sum of squared correlations outside blocks (log-scale)",
       x = "Sum of squared blocks", color = "Maximum block size") +
  geom_vline(xintercept = 0) +
  xlim(0, NA)

print(subset(res4, max_size == 4119), n = Inf)
