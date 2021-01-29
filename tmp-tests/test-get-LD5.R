# computed with snp_cor() (do not hesitate to be stringent on alpha)
corr <- readRDS(url("https://www.dropbox.com/s/65u96jf7y32j2mj/spMat.rds?raw=1"))

grid_param <- expand.grid(min_size = c(10, 20),
                          max_size = c(30, 40, 50),
                          lambda = c(0, 0.001, 0.01, 0.1))

# grid_param <- expand.grid(min_size = c(10, 15, 20, 25),
#                           max_size = c(30, 35, 40, 45, 50),
#                           lambda = 0)

THR_R2 <- 0.005

# debugonce(bigsnpr::snp_ldsplit)
system.time(
  res <- bigsnpr::snp_ldsplit(corr, grid_param, thr_r2 = THR_R2)
)
res[1:5]
#

#   min_size max_size lambda n_block    cost
# 1       10       30  0          18  8.74
# 2       20       30  0          15 10.2
# 3       10       40  0          15  3.14
# 4       20       40  0          14  3.14
# 5       10       50  0          13  0.0777
# 6       20       50  0          12  0.0777

# 7       10       30  0.001      18  8.74
# 8       20       30  0.001      15 10.2
# 9       10       40  0.001      15  3.14
# 10       20       40  0.001      14  3.14
# 11       10       50  0.001      13  0.0777
# 12       20       50  0.001      12  0.0777

# 13       10       30  0.01       20  8.79
# 14       20       30  0.01       15 10.3
# 15       10       40  0.01       19  3.39
# 16       20       40  0.01       14  3.14
# 17       10       50  0.01       17  0.240
# 18       20       50  0.01       13  0.127

best_res <- res[6, ]
all_ind <- head(best_res$all_last[[1]], -1)
block_num <- best_res$block_num[[1]]

library(dplyr)
library(ggplot2)

## Transform sparse representation into (i,j,x) triplets
corrT <- as(corr, "dgTMatrix")

ind_group <- split(seq_along(block_num), block_num)
ind_two <- lapply(seq(0, length(ind_group) - 2L), function(add) 1:2 + add)

STEP <- 20

all_plots <- lapply(ind_two, function(ind2) {

  ind <- unlist(ind_group[ind2], use.names = FALSE)

  df2 <- corrT[ind, ind] %>%
    { filter(data.frame(i = .@i + ind[1], j = .@j + ind[1], r2 = .@x^2), r2 >= THR_R2) }

  lims <- all_ind[ind2[1]] + 0.5

  range <- round(range(df2$i) / STEP)

  ggplot(df2, aes(x = i, y = j)) +
    geom_raster(aes(fill = r2)) +
    scale_fill_viridis_b(direction = -1, breaks = c(0.02, 0.1, 0.3, 0.8)) +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    geom_vline(xintercept = lims, linetype = 2) +
    geom_hline(yintercept = lims, linetype = 2) +
    scale_y_reverse(breaks = NULL) +
    scale_x_continuous(breaks = seq(range[1], range[2]) * STEP,
                       minor_breaks = (seq(range[1], range[2]) + 0.5) * STEP) +
    coord_equal() +
    theme(axis.text.y = element_blank(), strip.text.x = element_blank())
})

cowplot::plot_grid(
  plotlist = c(lapply(all_plots, function(p) p + theme(legend.position = "none")),
               list(cowplot::get_legend(all_plots[[1]]))),
  nrow = 3
)

# cowplot::plot_grid(
#   cowplot::plot_grid(
#     plotlist = c(all_plots[[1]] + theme(legend.position = "left"),
#                  lapply(all_plots[-1],
#                         function(p) p + theme(legend.position = "none"))),
#     nrow = 4
#   ),
#   cowplot::get_legend(all_plots[[1]]),
#   nrow = 1, rel_widths = c(7, 1)
# )

all_plots[[22]]
bigsnpr:::compute_cost(block_num, corr, THR_R2) # 4.79
block_num2 <- block_num; block_num2[401 - 0:20] <- 22
bigsnpr:::compute_cost(block_num2, corr, THR_R2) # 2.29
