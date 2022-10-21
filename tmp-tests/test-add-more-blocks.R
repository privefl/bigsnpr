corr <- readRDS("~/LD_with_blocks_chr22.rds")

info <- readRDS(url("https://figshare.com/ndownloader/files/36360900"))

block_id <- subset(info, chr == 22)$group_id - 595L
table(block_id)

ind <- 1:(271 + 50)
Matrix::image(corr[ind, ind]^2)

(ind <- which(block_id == 22))
corr2 <- corr[ind, ind]
corrT <- as(corr2 ** 2, "dgTMatrix")
df <- tibble::tibble(
  i = corrT@i + 1L,
  j = corrT@j + 1L,
  r2 = corrT@x
)
dim(df)
dim(df2 <- df[df$r2 > 0.01, ])

library(ggplot2)
ggplot(df2) +
  geom_tile(aes(i, j, color = r2, alpha = r2)) +
  coord_fixed() +
  theme_void() +
  scale_color_gradient(low = "#FFFFFF", high = "#000000", guide = "none") +
  scale_alpha(guide = "none")

library(ggplot2)
ggplot(subset(df, r2 > 0.01)) +
  geom_tile(aes(i, j, alpha = r2)) +
  coord_fixed() +
  theme_void() +
  geom_vline(xintercept = c(112, 220)) +
  scale_color_gradient(low = "#FFFFFF", high = "#000000", guide = "none") +
  scale_alpha(guide = "none")

compute_cost <- function(block_num, corr, thr_r2) {
  require("magrittr")
  corr %>%
    Matrix::tril() %>%
    as("dgTMatrix") %>%
    { .@x[block_num[.@i + 1L] != block_num[.@j + 1L]]^2 } %>%
    { .[. >= thr_r2] } %>%
    { c(sum(.), max(.)) }
}
first <- 106; second <- 239; compute_cost(rep(1:3, c(first, second - first, 300 - second)),
                                          corr2, thr_r2 = 0.01)

bigsnpr::snp_ldsplit(corr2, thr_r2 = 0.01, min_size = 10, max_size = 200,
                     max_r2 = 0.3, max_cost = 50)

bigsnpr::snp_ldsplit(corr2, thr_r2 = 0.01, min_size = 10, max_size = 100,
                     max_r2 = 0.3, max_cost = 1000) %>%
  print(n = Inf)

res <- .Last.value$all_size[[8]]
ggplot(subset(df, r2 > 0.01)) +
  geom_tile(aes(i, j, alpha = r2)) +
  coord_fixed() +
  theme_void() +
  geom_vline(xintercept = cumsum(head(res, -1)) + 0.5, linetype = 2) +
  geom_hline(yintercept = cumsum(head(res, -1)) + 0.5, linetype = 2) +
  scale_color_gradient(low = "#FFFFFF", high = "#000000", guide = "none") +
  scale_alpha(guide = "none")

bigsnpr::snp_ldsplit(corr2, thr_r2 = 0.02, min_size = 10, max_size = ncol(corr2) / 2,
                     max_r2 = 0.3, max_cost = Inf) %>%
  print(n = Inf)

res <- .Last.value$all_size[[3]]
ggplot(subset(df, r2 > 0.01)) +
  geom_tile(aes(i, j, alpha = r2^2)) +
  coord_fixed() +
  theme_void() +
  geom_vline(xintercept = cumsum(head(res, -1)) + 0.5, linetype = 2, color = "chartreuse") +
  geom_hline(yintercept = cumsum(head(res, -1)) + 0.5, linetype = 2, color = "chartreuse") +
  scale_color_gradient(low = "#FFFFFF", high = "#000000", guide = "none") +
  scale_alpha(guide = "none")


corr2 <- corr[ind, ind]
corr3 <- Matrix::nearPD(0.95 * corr2 + Matrix::Diagonal(n = ncol(corr2), x = 0.05))
corr3$mat[1:5, 1:5]

library(bigstatsr)
df2 <- expand.grid(i = rows_along(corr2), j = cols_along(corr2))
df2$r2 <- as.vector(as.matrix(corr3$mat))
ggplot(df2) +
  geom_tile(aes(i, j, color = r2, alpha = r2)) +
  coord_fixed() +
  theme_void() +
  scale_color_gradient(low = "#FFFFFF", high = "#000000", guide = "none") +
  scale_alpha(guide = "none")

plot(corr3$mat[, 1])
