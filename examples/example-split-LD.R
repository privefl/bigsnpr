\dontrun{

  corr <- readRDS(url("https://www.dropbox.com/s/65u96jf7y32j2mj/spMat.rds?raw=1"))

  # adjust `THR_R2` depending on sample size used to compute corr
  # use e.g. 0.05 for small sample sizes, and 0.01 for large sample sizes
  THR_R2 <- 0.02
  m <- ncol(corr)
  (SEQ <- round(seq_log(m / 30, m / 5, length.out = 10)))
  # replace `min_size` by e.g. 100 for larger data
  (res <- snp_ldsplit(corr, thr_r2 = THR_R2, min_size = 10, max_size = SEQ))

  # add the variant block IDs corresponding to each split
  res$block_num <- lapply(res$all_size, function(.) rep(seq_along(.), .))

  library(ggplot2)
  # trade-off cost / number of blocks
  qplot(n_block, cost, color = factor(max_size, SEQ), data = res) +
    theme_bw(14) +
    scale_y_log10() +
    theme(legend.position = "top") +
    labs(x = "Number of blocks", color = "Maximum block size",
         y = "Sum of squared correlations outside blocks")

  # trade-off cost / number of non-zero values
  qplot(perc_kept, cost, color = factor(max_size, SEQ), data = res) +
    theme_bw(14) +
    # scale_y_log10() +
    theme(legend.position = "top") +
    labs(x = "Percentage of non-zero values kept", color = "Maximum block size",
         y = "Sum of squared correlations outside blocks")

  # trade-off cost / sum of squared sizes
  qplot(cost2, cost, color = factor(max_size, SEQ), data = res) +
    theme_bw(14) +
    scale_y_log10() +
    geom_vline(xintercept = 0)+
    theme(legend.position = "top") +
    labs(x = "Sum of squared blocks", color = "Maximum block size",
         y = "Sum of squared correlations outside blocks")


  ## Pick one solution and visualize blocks
  library(dplyr)
  all_ind <- res %>%
    arrange(cost2 * sqrt(5 + cost)) %>%
    print() %>%
    slice(1) %>%
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
