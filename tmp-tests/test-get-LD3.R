library(bigsnpr)
bedfile <- download_1000G("tmp-data")
obj.bed <- bed(bedfile)
fam2 <- bigreadr::fread2(sub_bed(bedfile, ".fam2"))
ind.eur <- which(fam2$`Super Population` == "EUR")
map <- bigreadr::fread2(sub_bed(bedfile, ".bim"))
chr <- 12
ind.chr <- which(map$V1 == chr)
keep_mac <- bed_MAF(obj.bed, ind.row = ind.eur, ind.col = ind.chr)$mac > 10
rds <- snp_readBed2(bedfile, tempfile(),
                    ind.row = ind.eur, ind.col = ind.chr[keep_mac])
bigsnp <- snp_attach(rds)
G <- bigsnp$genotypes
POS <- bigsnp$map$physical.pos
POS2 <- snp_asGeneticPos(rep(chr, length(POS)), POS, dir = "tmp-data")

corr <- snp_cor(G, infos.pos = POS2, size = 3 / 1000, ncores = 6)
corr[1:5, 1:5]

grid_param <- expand.grid(min_size = c(500, 1000),
                          max_size = c(5000, 10e3),
                          lambda = c(0, 1e-3))

grid_param <- expand.grid(min_size = c(1000),
                          max_size = c(10e3, 14e3, 18e3),
                          lambda = c(0, 1e-2, 1e-3))


system.time(
  res <- snp_ldsplit(corr, grid_param, thr_r2 = 0.05)
)
# 40 sec for chr 22 for first grid
# 5 min for chr 6 for second grid

res
# cost was 1007 for chr 22 -> now is 41

#   min_size max_size lambda n_block  cost
# 1     1000    10000  0          13 7283.
# 2     1000    14000  0           9  613.
# 3     1000    18000  0           7  327.
# 4     1000    10000  0.01       15 7360.
# 5     1000    14000  0.01       12  757.
# 6     1000    18000  0.01       11  569.
# 7     1000    10000  0.001      13 7283.
# 8     1000    14000  0.001      10  613.
# 9     1000    18000  0.001       7  327.

# min_size max_size lambda n_block    cost
# 1     1000    10000  0          13 1398.
# 2     1000    14000  0          10   18.0
# 3     1000    18000  0           8    5.55
# 4     1000    10000  0.01       48 1570.
# 5     1000    14000  0.01       43  181.
# 6     1000    18000  0.01       41  179.
# 7     1000    10000  0.001      27 1413.
# 8     1000    14000  0.001      24   34.6
# 9     1000    18000  0.001      22   25.6

all_ind <- head(res$all_last[[8]], -1)

# PLOT
{
  corrT <- as(corr ** 2, "dgTMatrix")
  upper <- which((corrT@i <= corrT@j) & (abs(corrT@x) > 0.05))
  # upper <- sample(upper, min(300e3, length(upper)))
  POS3 <- seq_len(ncol(corr))
  df0 <- tibble::tibble(
    i = POS3[corrT@i[upper] + 1L],
    j = POS3[corrT@j[upper] + 1L],
    r2 = corrT@x[upper],
    y = (j - i) / 2,
    z = i + y
  )
  # hist(df0$y)
  # dim(df <- dplyr::slice_max(df0, y, n = 200e3))
  df <- dplyr::slice_max(df0, y, n = 300e3)

  K <- ceiling(max(POS2) / 30)
  breaks <- quantile(range(df$z), probs = 0:K / K)
  breaks[1] <- breaks[1] - 1
}


library(ggplot2)
library(dplyr)
ggplot(mutate(df, cut = cut(z, breaks))) +
  geom_point(aes(z, y, color = r2), size = rel(0.5)) +
  # coord_fixed() +
  scale_color_gradientn(colours = rev(colorRamps::matlab.like2(100)), guide = "none") +
  theme_minimal() +
  scale_x_continuous(expand = c(0.02, 0.02)) +
  geom_vline(aes(xintercept = z), linetype = 3,
             data = tibble(z = (POS3[all_ind] + POS3[all_ind + 1L]) / 2,
                           cut = cut(z, breaks))) +
  scale_y_continuous(breaks = 0, limits = c(0, NA)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  labs(x = "Position", y = NULL) +
  # theme(aspect.ratio = 1) +
  scale_alpha(guide = "none") + #scale_color_continuous(guide = "none") +
  facet_wrap(~ cut, scales = "free_x", ncol = 1) +
  theme(strip.background = element_blank(), strip.text.x = element_blank())

# POS[all_ind[13:14]] # HLA: 25-34
# diff(POS2[all_ind[13:14] + 1:0]) # 3.96
