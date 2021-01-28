\dontrun{

  corr <- readRDS(url("https://www.dropbox.com/s/65u96jf7y32j2mj/spMat.rds?raw=1"))

  grid_param <- expand.grid(min_size = c(10, 20),
                            max_size = c(30, 40, 50),
                            lambda = c(0, 1e-3, 1e-2))
  THR_R2 <- 0.005

  (res <- snp_ldsplit(corr, grid_param, thr_r2 = THR_R2))

  all_ind <- head(res$all_last[[5]], -1)

  ## Transform sparse representation into (i,j,x) triplets
  corrT <- as(corr, "dgTMatrix")
  upper <- (corrT@i <= corrT@j & corrT@x^2 >= THR_R2)
  df <- data.frame(
    i = corrT@i[upper] + 1L,
    j = corrT@j[upper] + 1L,
    r2 = corrT@x[upper]^2
  )
  df$y <- (df$j - df$i) / 2

  library(ggplot2)
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
