library(bigsnpr)
bedfile <- download_1000G("tmp-data")
obj.bed <- bed(bedfile)
fam2 <- bigreadr::fread2(sub_bed(bedfile, ".fam2"))
ind.eur <- which(fam2$`Super Population` == "EUR")
map <- bigreadr::fread2(sub_bed(bedfile, ".bim"))
chr <- 15
ind.chr <- which(map$V1 == chr)
keep_mac <- bed_MAF(obj.bed, ind.row = ind.eur, ind.col = ind.chr)$mac > 10
rds <- snp_readBed2(bedfile, tempfile(),
                    ind.row = ind.eur, ind.col = ind.chr[keep_mac])
bigsnp <- snp_attach(rds)
G <- bigsnp$genotypes
POS <- bigsnp$map$physical.pos
POS2 <- snp_asGeneticPos(rep(chr, length(POS)), POS, dir = "tmp-data")

corr <- snp_cor(G, infos.pos = POS2, size = 4 / 1000, ncores = 6)
corr[1:5, 1:5]

library(dplyr)
corrT <- as(corr, "dgTMatrix")
upper <- (corrT@i <= corrT@j)
m <- ncol(corr)
all_ind <- data.frame(
  i = factor(1:m)[corrT@i[upper] + 1L],
  j = corrT@j[upper] + 1L,
  r2 = corrT@x[upper] ** 2
) %>%
  filter(r2 > 0.1) %>%
  group_by(i) %>%
  summarize(max_j = max(j)) %>%
  mutate(cum_max_j = cummax(max_j)) %>%
  filter(i == cum_max_j) %>%
  pull(i) %>%
  print()

MIN_M <- 1000
MAX_M <- 5000

# lower_ld <- Matrix::colSums(Matrix::tril(corr, -1) ** 2)
# lower_ld_smooth <- bigutilsr::rollmean(lower_ld, 10)
# plot(lower_ld_smooth, pch = 20, cex = 0.5)

library(Matrix)

m <- ncol(corr)

Rcpp::sourceCpp('src/split-LD.cpp')

library(magrittr)
E <- corr %>%
  Matrix::tril() %>%
  { get_L(.@p, .@i, .@x, thr_r2 = 0.05) } %>%  # res
  { Matrix::sparseMatrix(i = .[[1]], j = .[[2]], x = .[[3]], dims = c(m, m),
                         triangular = FALSE, index1 = FALSE) } %>%  # L
  get_E(min_row = pmin(min_row(corr@p, corr@i), pmax(1:m - MAX_M, 0))) %>%
  subset((j - i + 1) > MIN_M) %>%  # res2
  { split(.[c("j", "x")], factor(1:m)[.$i]) }

str(E)

LAMBDA <- 0 #1e-3

best_ind <- rep(NA, m)
C <- rep(NA, m)
C[m - MIN_M:MAX_M + 1] <- 0
for (i in m:1) {
  e <- E[[i]]
  cost <- e$x + C[e$j + 1L] + LAMBDA * (e$j - i)
  ind.min <- which.min(cost)
  if (length(ind.min) > 0) {
    best_ind[i] <- e$j[ind.min]
    C[i] <- e$x[ind.min] + C[e$j[ind.min] + 1L]
  }
}
plot(C)

# reconstruct path
all_ind <- list()
j <- 0
repeat {
  j <- best_ind[j + 1L]
  print(j)
  if (is.na(j)) break
  all_ind[[length(all_ind) + 1L]] <- j
}
all_ind <- unlist(all_ind)
# all_ind <- all_ind[(m - all_ind) > MIN_M]


# PLOT
{
  corrT <- as(corr ** 2, "dgTMatrix")
  upper <- which((corrT@i <= corrT@j) & (abs(corrT@x) > 0.05))
  upper <- sample(upper, min(500e3, length(upper)))
  POS3 <- seq_len(m)
  df0 <- tibble::tibble(
    i = POS3[corrT@i[upper] + 1L],
    j = POS3[corrT@j[upper] + 1L],
    r2 = corrT@x[upper],
    y = (j - i) / 2,
    z = i + y
  )
  hist(df0$y)
  dim(df <- dplyr::slice_max(df0, y, n = 200e3))

  K <- ceiling(max(POS2) / 25)
  breaks <- quantile(range(df$z), probs = 0:K / K)
  breaks[1] <- breaks[1] - 1
}


library(ggplot2)
library(dplyr)
ggplot(mutate(df, cut = cut(z, breaks))) +
  geom_point(aes(z, y, color = r2, alpha = r2), size = rel(0.5)) +
  # coord_fixed() +
  scale_color_gradientn(colours = rev(colorRamps::matlab.like2(100)), guide = "none") +
  theme_minimal() +
  scale_x_continuous(expand = c(0.02, 0.02)) +
  geom_vline(aes(xintercept = z), linetype = 3,
             data = tibble(z = POS3[all_ind], cut = cut(z, breaks))) +
  scale_y_continuous(breaks = 0, limits = c(0, NA)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  labs(x = "Position", y = NULL) +
  # theme(aspect.ratio = 1) +
  scale_alpha(guide = "none") + #scale_color_continuous(guide = "none") +
  facet_wrap(~ cut, scales = "free_x", ncol = 1) +
  theme(strip.background = element_blank(), strip.text.x = element_blank())

# POS[all_ind[13:14]] # HLA: 25-34
# diff(POS2[all_ind[13:14] + 1:0]) # 3.96
c(C[1], length(all_ind), MIN_M, MAX_M, LAMBDA)
# 24 / 18 / 1000 / 5000 / 1e-3
# 19 / 11 / 1000 / 5000 / 0
