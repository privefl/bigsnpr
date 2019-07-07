bedfile <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
# bedfile <- "../Dubois2010_data/FinnuncorrNLITUK1UK3hap300_QC_norel.bed"
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

system.time(bed_autoSVD(bedfile, ncores = nb_cores()))
