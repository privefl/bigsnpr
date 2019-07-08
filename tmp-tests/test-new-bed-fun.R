library(bigsnpr)
bedfile <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
bedfile <- "../Dubois2010_data/FinnuncorrNLITUK1UK3hap300_QC_norel.bed"
assert_exist <- bigsnpr:::assert_exist
NAMES.MAP <- bigsnpr:::NAMES.MAP
NAMES.FAM <- bigsnpr:::NAMES.FAM
obj.bed <- bed(bedfile)
obj.bed$bedfile
obj.bed$famfile
obj.bed$nrow
obj.bed$ncol
ind.row <- rows_along(obj.bed)
ind.col <- cols_along(obj.bed)

system.time(ind.keep <- bed_clumping(bedfile, ncores = nb_cores(),
                                     exclude = bed_indLRLDR(bedfile))) # 55-60 sec

stats <- bigsnpr:::bed_stats(obj.bed, ind.row, ind.col)
af <- stats$sum / (2 * stats$nb_nona_col)
center <- 2 * af
scale <- sqrt(2 * af * (1 - af))

bigsnpr:::cpMatVec4(obj.bed, ind.row, ind.col, center, scale,
                    rnorm(obj.bed$nrow))

system.time(svd <- bed_randomSVD(bedfile, ind.col = ind.keep))  # 230 sec -> 500
plot(svd)
plot(svd, type = "scores")
plot(svd, type = "scores", scores = 3:4)
plot(svd, type = "scores", scores = 5:6)
plot(svd, type = "loadings", loadings = 1:6, coeff = 0.5)

cor(stats$nb_nona_row, svd$u)

system.time(svd2 <- bed_randomSVD(bedfile, ind.col = ind.keep,
                                  ncores = nb_cores()))  # 184 sec -> 356
all.equal(svd2, svd)

system.time(svd3 <- bed_autoSVD(bedfile, ncores = nb_cores(), thr.r2 = 0.2))
# 2091 sec -> 35 min
attr(svd3, "lrldr")
plot(svd3$u[, 1:3], svd2$u[, 1:3], pch = 20); abline(0, 1, col = "red")
plot(svd3)
plot(svd3, type = "loadings", loadings = 1:6, coeff = 0.5)
plot(svd3, type = "scores")
plot(svd3, type = "scores", scores = 3:4)
obj.bed <- bed(bedfile)
pop <- obj.bed$fam$family.ID
library(ggplot2)
plot(svd3, type = "scores", scores = 3:4) +
  aes(color = pop)

CHR <- obj.bed$map$chromosome
POS <- obj.bed$map$physical.pos
gwas <- bed_pcadapt(bedfile, U.row = svd3$u, ncores = nb_cores())
snp_manhattan(gwas, CHR, POS, npoints = 20e3) +
  geom_hline(yintercept = -log10(5e-8), color = "red")
snp_qq(gwas) + xlim(1, NA)

snp_manhattan(gwas, CHR, POS, npoints = 20e3,
              ind.highlight = snp_indLRLDR(CHR, POS, LD.regions = rbind(attr(svd3, "lrldr"), LD.wiki34[1:3]))) +
  geom_hline(yintercept = -log10(5e-8), color = "red")

snp_manhattan(gwas, CHR, POS, npoints = 20e3,
              ind.highlight = attr(svd3, "subset")) +
  geom_hline(yintercept = -log10(5e-8), color = "red")

ind.keep2 <- setdiff(attr(svd3, "subset"), which(-predict(gwas) > 5))
system.time(
  svd4 <- bed_randomSVD(bedfile, ind.col = ind.keep2, ncores = nb_cores())
)

gwas2 <- bed_pcadapt(bedfile, U.row = svd4$u, ncores = nb_cores())
snp_manhattan(gwas2, CHR, POS, npoints = 20e3,
              ind.highlight = ind.keep2) +
  geom_hline(yintercept = -log10(5e-8), color = "blue")

plot(svd4, type = "scores")
plot(svd4, type = "scores", scores = 3:4)
plot(svd4, type = "scores", scores = 5:6)
plot(svd4, type = "scores", scores = 7:8)
plot(svd4, type = "scores", scores = 9:10)

devtools::source_gist("42b41d771bbeae63245b8304ef283c70", filename = "get-genes.R")
rsid <- obj.bed$map$marker.ID[which(predict(gwas2, log10 = FALSE) < 5e-8)]
genes <- snp_gene(rsid)
genes

dist <- robust::covRob(svd$u, estim = "pairwiseGK")$dist
plot(svd4, type = "scores") +
  coord_fixed() +
  aes(color = log(dist)) +
  scale_colour_viridis_c()
plot(svd4, type = "scores", scores = 3:4) +
  coord_fixed() +
  aes(color = log(dist)) +
  scale_colour_viridis_c()
plot(svd4, type = "scores", scores = 5:6) +
  coord_fixed() +
  aes(color = log(dist)) +
  scale_colour_viridis_c()
plot(svd4, type = "scores", scores = 7:8) +
  coord_fixed() +
  aes(color = log(dist)) +
  scale_colour_viridis_c()
plot(svd4, type = "scores", scores = 9:10) +
  coord_fixed() +
  aes(color = log(dist)) +
  scale_colour_viridis_c()

ind.row <- which(log(dist) < 4 & obj.bed$fam$affection == 1)

system.time(
  svd5 <- bed_autoSVD(bedfile, ind.row = ind.row, ind.col = ind.keep2, ncores = nb_cores())
)

plot(svd5, type = "scores")
plot(svd5, type = "scores", scores = 3:4)
plot(svd5, type = "scores", scores = 5:6)

dist2 <- robust::covRob(svd5$u, estim = "pairwiseGK")$dist
plot(svd5, type = "scores") +
  coord_fixed() +
  aes(color = log(dist2), alpha = -log(dist2)) +
  scale_colour_viridis_c()
plot(svd5, type = "scores", scores = 3:4) +
  coord_fixed() +
  aes(color = log(dist2)) +
  scale_colour_viridis_c()

ind.row2 <- ind.row[log(dist2) < 3]
system.time(
  svd6 <- bed_autoSVD(bedfile, ind.row = ind.row2,
                      ind.col = ind.keep2, ncores = nb_cores())
)
plot(svd6, type = "scores")
plot(svd6, type = "scores", scores = 3:4)
