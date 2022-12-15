corr <- readRDS("~/LD_with_blocks_chr22.rds")

library(dplyr)
res0 <- bigsnpr::snp_ldsplit(corr^2, thr_r2 = 0.01,
                             min_size = 8, max_size = 100,
                             max_K = ceiling(ncol(corr) / 30),
                             max_r2 = 1, max_cost = Inf) %>%
  print(n = Inf)

library(ggplot2)
qplot(perc_kept, cost, color = factor(max_size), data = res0) +
  theme_bw(14) +
  scale_y_log10() +
  theme(legend.position = "top") +
  labs(x = "Percentage of non-zero values kept", color = "Maximum block size",
       y = "Sum of squared correlations outside blocks")


res <- res0 %>%
  # filter(cost < (1.1 * min(cost))) %>%
  # filter(perc_kept < 0.025) %>%
  # slice_min(cost) %>%
  slice_max(n_block) %>%
  pull(all_last) %>%
  .[[1]] %>%
  print()

corr_to_df <- function(corr, ind, thr = 0.02) {
  corr[ind, ind] %>%
    Matrix::summary() %>%
    as.data.frame() %>%
    dplyr::transmute(i = ind[i], j = ind[j], r2 = x^2) %>%
    dplyr::filter(r2 > thr)
}

ind <- 300:500
ind <- 500:1000 + 1000
ggplot(corr_to_df(corr2, ind)) +
  geom_tile(aes(i, j, alpha = r2)) +
  coord_fixed() +
  theme_void() +
  geom_vline(xintercept = res[res %in% ind] + 0.5, linetype = 3, color = "red") +
  geom_hline(yintercept = res[res %in% ind] + 0.5, linetype = 3, color = "red") +
  scale_color_gradient(low = "#FFFFFF", high = "#000000", guide = "none") +
  scale_alpha(guide = "none")

all_size <- diff(c(0, res))
block_num <- rep(seq_along(all_size), all_size)
ind_block <- split(seq_along(block_num), block_num)

corr2 <- as(corr, "dgCMatrix")
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
      stat <- sqrt(Matrix::mean(corr2[ind_i, ind_j]^4))
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

library(foreach)
ind_keep <- foreach(j = 1:K, .combine = "rbind") %do% {
  lim <- range(which(keep[, j]))
  ind_j <- ind_block[[j]]
  ind_i <- seq(head(ind_block[[lim[1]]], 1), tail(ind_block[[lim[2]]], 1))
  expand.grid(ind_i, ind_j)
}

corr3 <- Matrix::sparseMatrix(
  i = ind_keep[[1]], j = ind_keep[[2]], x = corr2[as.matrix(ind_keep)],
  dims = dim(corr2), index1 = TRUE, symmetric = FALSE)
class(corr3)
Matrix::isSymmetric(corr3)

corr4 <- as(Matrix::drop0(corr3), "symmetricMatrix")
# c(length(corr2@x), length(corr3@x), length(corr4@x))
# diff <- corr3 - corr2
# hist(diff@x)
# summary(diff@x)
length(corr4@x) / length(corr@x) # 13% // 9%
# Matrix::image(corr4^2)

ind2 <- 500:1000 + 1000
# ind2 <- ind
ggplot(corr_to_df(corr2, ind2), aes(i, j, alpha = r2)) +
  geom_tile(data = corr_to_df(corr3 != 0, ind2), color = "skyblue", alpha = 0.001) +
  geom_tile() +
  coord_fixed() +
  theme_void() +
  geom_vline(xintercept = res[res %in% ind2] + 0.5, linetype = 3, color = "red") +
  geom_hline(yintercept = res[res %in% ind2] + 0.5, linetype = 3, color = "red") +
  scale_color_gradient(low = "#FFFFFF", high = "#000000", guide = "none") +
  scale_alpha(guide = "none")

RSpectra::eigs(corr0, k = 2)$values
# 106.13132  83.42977
RSpectra::eigs(corr0, k = 2, sigma = -0.05)$values
# -0.007954082 -0.009787068

RSpectra::eigs(corr4, k = 2)$values
# 104.45772  76.84316
RSpectra::eigs(corr4, k = 2, sigma = -10)$values
# -1.234222 -1.261020


prec <- solve(corr0[1:5, 1:5])
(pcor <- -cov2cor(prec))
diag(pcor) <- 1
corr0[1:5, 1:5]
(pcor2 <- corpcor::cor2pcor(corr0[1:5, 1:5]))
all.equal(pcor, pcor2)
