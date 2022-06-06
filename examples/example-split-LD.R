\dontrun{

  corr <- readRDS(url("https://www.dropbox.com/s/65u96jf7y32j2mj/spMat.rds?raw=1"))

  THR_R2 <- 0.01

  # You can provide multiple max_size to try
  # For real data, use e.g. round(seq_log(1000 + ncol(corr) / 30,
  #                                       5000 + ncol(corr) / 10,
  #                                       length.out = 6))
  SEQ <- c(50, 80)
  (res <- snp_ldsplit(corr, thr_r2 = THR_R2, min_size = 10, max_size = SEQ,
                      max_K = 50, max_r2 = 0.5, max_cost = Inf))

  library(ggplot2)
  # trade-off cost / number of blocks
  qplot(n_block, cost, color = factor(max_size, SEQ), data = res) +
    theme_bw(14) +
    scale_y_log10() +
    theme(legend.position = c(0.5, 0.82)) +
    labs(x = "Number of blocks", color = "Maximum block size",
         y = "Sum of squared correlations outside blocks")

  # trade-off cost / number of non-zero values
  qplot(perc_kept, cost, color = factor(max_size, SEQ), data = res) +
    theme_bw(14) +
    # scale_y_log10() +
    theme(legend.position = c(0.6, 0.82)) +
    labs(x = "Percentage of non-zero values kept", color = "Maximum block size",
         y = "Sum of squared correlations outside blocks")

  # trade-off cost / sum of squared sizes
  res$cost2 <- sapply(res$all_size, function(sizes) sum(sizes^2))
  qplot(cost2, cost, color = factor(max_size, SEQ), data = res) +
    theme_bw(14) +
    scale_y_log10() +
    geom_vline(xintercept = 0)+
    theme(legend.position = c(0.7, 0.82)) +
    labs(x = "Sum of squared blocks", color = "Maximum block size",
         y = "Sum of squared correlations outside blocks")


  ## Pick one solution and visualize blocks
  library(dplyr)
  all_ind <- res %>%
    filter(cost < 1) %>%
    arrange(perc_kept) %>%
    slice(1) %>%
    print() %>%
    pull(all_last)

  ## Transform sparse representation into (i,j,x) triplets
  corrT <- as(corr, "dgTMatrix")
  upper <- (corrT@i <= corrT@j & corrT@x^2 >= THR_R2)
  df <- data.frame(
    i = corrT@i[upper] + 1L,
    j = corrT@j[upper] + 1L,
    r2 = corrT@x[upper]^2
  )
  df$y <- (df$j - df$i) / 2

  ggplot(df) +
    geom_point(aes(i + y, y, alpha = r2)) +
    theme_minimal() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          strip.background = element_blank(), strip.text.x = element_blank()) +
    scale_alpha_continuous(range = 0:1) +
    scale_x_continuous(expand = c(0.02, 0.02), minor_breaks = NULL,
                       breaks = head(all_ind[[1]], -1) + 0.5) +
    facet_wrap(~ cut(i + y, 4), scales = "free", ncol = 1) +
    labs(x = "Position", y = NULL)
}
