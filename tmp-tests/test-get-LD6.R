library(dplyr)
library(ggplot2)

corr <- readRDS(runonce::download_file(
  "https://ndownloader.figshare.com/files/24928052", dir = "tmp-data"))

## Transform sparse representation into (i,j,x) triplets
corrT <- as(corr, "dgTMatrix")

grid_param <- expand.grid(min_size = c(200),
                          max_size = c(2000, 3000, 4000),
                          lambda = c(0, 0.001, 0.003, 0.01, 0.03))
THR_R2 <- 0.005

# debugonce(bigsnpr::snp_ldsplit)
system.time(
  res <- bigsnpr::snp_ldsplit(corr, grid_param, thr_r2 = THR_R2)
)
print(res[1:5], n = 30)

qplot(n_block, cost, data = res) + theme_bw(16)

best_res <- res[13, ]
all_ind <- head(best_res$all_last[[1]], -1)
block_num <- best_res$block_num[[1]]

ind_group <- split(seq_along(block_num), block_num)
ind_two <- lapply(seq(0, length(ind_group) - 2L), function(add) 1:2 + add)

STEP <- 500

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

all_plots[[10]]  # why not more even blocks?
bigsnpr:::compute_cost(block_num, corr, THR_R2) # 4.79
block_num2 <- block_num; block_num2[401 - 0:20] <- 22
bigsnpr:::compute_cost(block_num2, corr, THR_R2) # 2.29
