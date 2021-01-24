library(bigsnpr)
bedfile <- download_1000G("tmp-data")
obj.bed <- bed(bedfile)
fam2 <- bigreadr::fread2(sub_bed(bedfile, ".fam2"))
ind.eur <- which(fam2$`Super Population` == "EUR")
map <- bigreadr::fread2(sub_bed(bedfile, ".bim"))
chr <- 6
ind.chr <- which(map$V1 == chr)
keep_mac <- bed_MAF(obj.bed, ind.row = ind.eur, ind.col = ind.chr)$mac > 10
rds <- snp_readBed2(bedfile, tempfile(),
                    ind.row = ind.eur, ind.col = ind.chr[keep_mac])
bigsnp <- snp_attach(rds)
G <- bigsnp$genotypes
POS <- bigsnp$map$physical.pos
POS2 <- snp_asGeneticPos(rep(chr, length(POS)), POS, dir = "tmp-data")

corr <- snp_cor(G, infos.pos = POS2, size = 3 / 1000, ncores = 4)
corr[1:5, 1:5]

# lower_ld <- Matrix::colSums(Matrix::tril(corr, -1) ** 2)
# lower_ld_smooth <- bigutilsr::rollmean(lower_ld, 10)
# plot(lower_ld_smooth, pch = 20, cex = 0.5)

library(Matrix)
corr2 <- tril(corr, -1)
x <- corr2@x
ind <- which(x^2 < 0.05)
length(ind) / length(corr2@x)
x[ind] <- 0
Rcpp::sourceCpp('tmp-tests/test-getL.cpp')
res <- get_L(corr2@p, corr2@i, x)
L <- sparseMatrix(i = res[[1]], j = res[[2]], x = res[[3]], dims = dim(corr2),
                  triangular = FALSE, index1 = FALSE)
L[1:5, 1:8]

res2 <- subset(get_E(L@p, L@i, L@x), (j - i) > 200)
library(dplyr)
# Could return this directly from C
m <- ncol(corr)
E <- split(res2[c("j", "x")], factor(1:m)[res2$i])
str(E)

best_ind <- rep(NA, m)
C <- rep(NA, m)
C[i <- m] <- 0
repeat {
  print(i <- i - 1L)
  e <- E[[i]]
  if (nrow(e) > 0) {
    cost <- e$x + C[e$j + 1L]
    ind.min <- which.min(cost)
    best_ind[i] <- e$j[ind.min]
    C[i] <- cost[ind.min]
  } else {
    C[i] <- 0
  }
  if (i == 1L) break
}

all_ind <- list(i <- 1)
repeat {
  (i <- best_ind[i] + 1L)
  if (is.na(i)) break
  all_ind[[length(all_ind) + 1L]] <- i
}
all_ind <- unlist(all_ind)


# PLOT
corrT <- as(corr ** 2, "dgTMatrix")
upper <- which((corrT@i <= corrT@j) & (abs(corrT@x) > 0.05))
upper <- sample(upper, 500e3)
POS3 <- POS2
df0 <- tibble::tibble(
  i = POS3[corrT@i[upper] + 1L],
  j = POS3[corrT@j[upper] + 1L],
  r2 = corrT@x[upper],
  y = (j - i) / 2,
  z = i + y
)
hist(df0$y)
dim(df <- dplyr::slice_max(df0, y, n = 100e3))

K <- 6
breaks <- quantile(range(df$z), probs = 0:K / K)
breaks[1] <- breaks[1] - 1


library(ggplot2)
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

