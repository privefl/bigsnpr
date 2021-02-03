\dontrun{

  corr <- readRDS(url("https://www.dropbox.com/s/65u96jf7y32j2mj/spMat.rds?raw=1"))

  THR_R2 <- 0.01

  (res <- snp_ldsplit(corr, thr_r2 = THR_R2, min_size = 10, max_size = 50, max_K = 50))

  library(ggplot2)
  qplot(n_block, cost, data = res) + theme_bw(16) + scale_y_log10()

  all_ind <- head(res$all_last[[6]], -1)

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
    geom_point(aes(i + y, y, color = r2), size = rel(0.5)) +
    coord_fixed() +
    scale_color_gradientn(colours = rev(colorRamps::matlab.like2(100))) +
    theme_minimal() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    geom_vline(xintercept = all_ind + 0.5, linetype = 3) +
    labs(x = "Position", y = NULL) +
    scale_alpha(guide = 'none')
}
