corr <- readRDS("~/LD_with_blocks_chr9.rds")

info <- readRDS(runonce::download_file(
  "https://figshare.com/ndownloader/files/37802721",
  dir = "tmp-data", fname = "map_hm3_plus.rds"))

info_chr <- subset(info, chr == 9)

block_id <- info_chr$block_id - 154L
table(block_id)

(ind <- which(block_id == 2))
corr2 <- corr[ind, ind]
corrT <- as(corr2 ** 2, "dgTMatrix")
df <- tibble::tibble(
  i = corrT@i + 1L,
  j = corrT@j + 1L,
  r2 = corrT@x
)
dim(df)
nrow(df) / length(corrT)
dim(df2 <- df[df$r2 > 0.1, ])

library(ggplot2)
ggplot(df2) +
  geom_tile(aes(i, j, alpha = r2)) +
  coord_fixed() +
  theme_void() +
  # scale_color_gradient(low = "#FFFFFF", high = "#000000", guide = "none") +
  scale_alpha(guide = "none")


library(dplyr)

res0 <- bigsnpr::snp_ldsplit(corr2, thr_r2 = 0.8,
                             min_size = 8, max_size = 200, max_K = 1000,
                             max_r2 = 0.95, max_cost = Inf) %>%
  print(n = Inf)

qplot(perc_kept, cost, color = factor(max_size), data = res0) +
  theme_bw(14) +
  # scale_y_log10() +
  theme(legend.position = "top") +
  labs(x = "Percentage of non-zero values kept", color = "Maximum block size",
       y = "Sum of squared correlations outside blocks")


res <- res0 %>%
  filter(perc_kept < 0.025) %>%
  slice_min(cost) %>%
  pull(all_last) %>%
  .[[1]] %>%
  print()

ggplot(df2) +
  geom_tile(aes(i, j, alpha = r2)) +
  coord_fixed() +
  theme_void() +
  geom_vline(xintercept = head(res, -1) + 0.5, linetype = 3, color = "red") +
  geom_hline(yintercept = head(res, -1) + 0.5, linetype = 3, color = "red") +
  scale_color_gradient(low = "#FFFFFF", high = "#000000", guide = "none") +
  scale_alpha(guide = "none") +
  xlim(NA, 500) + ylim(NA, 500)

all_size <- diff(c(0, res))
block_num <- rep(seq_along(all_size), all_size)
ind_block <- split(seq_along(block_num), block_num)

corr3 <- as(corr2, "dgCMatrix")
K <- length(ind_block)
keep <- matrix(FALSE, K, K); diag(keep) <- TRUE
for (i in 1:K) {
  stop_next <- FALSE
  for (j in 1:K) {
    if (j <= i) {
      next
    } else {
      ind_i <- ind_block[[i]]
      ind_j <- ind_block[[j]]
      stat <- sqrt(Matrix::mean(corr3[ind_i, ind_j]^4))
      print(round(c(i, j, 100 * stat), 1))
      if (stat < 0.015) {
        if (stop_next) {
          break
        } else {
          stop_next <- TRUE
        }
      } else {
        keep[i:j, i:j] <- TRUE
        stop_next <- FALSE
      }
    }
  }
}


for (j in 1:K) {
  print(j)
  lim <- range(which(keep[, j]))
  ind_j <- ind_block[[j]]
  ind_i <- seq(head(ind_block[[lim[1]]], 1), tail(ind_block[[lim[2]]], 1))
  corr3[-ind_i, ind_j] <- 0
}

corr4 <- as(Matrix::drop0(corr3), "symmetricMatrix")
# c(length(corr2@x), length(corr3@x), length(corr4@x))
# diff <- corr3 - corr2
# hist(diff@x)
# summary(diff@x)
length(corr4@x) / length(corr2@x) # 13% // 9%
# Matrix::image(corr4^2)

ind2 <- 500:1000 + 1000
df3 <- Matrix::which(corr4 != 0, arr.ind = TRUE)
ggplot(filter(df2, i %in% ind2, j %in% ind2)) +
  geom_tile(aes(row, col),
            data = filter(as.data.frame(df3), row %in% ind2, col %in% ind2),
            color = "skyblue", alpha = 0.01) +
  geom_tile(aes(i, j, alpha = r2)) +
  coord_fixed() +
  theme_void() +
  geom_vline(xintercept = res[res %in% ind2] + 0.5, linetype = 3, color = "red") +
  geom_hline(yintercept = res[res %in% ind2] + 0.5, linetype = 3, color = "red") +
  scale_color_gradient(low = "#FFFFFF", high = "#000000", guide = "none") +
  # xlim(NA, 600) + ylim(NA, 600) +
  scale_alpha(guide = "none")
