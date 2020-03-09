library(bigsnpr)

celiac <- snp_attach("../Dubois2010_data/FinnuncorrNLITUK1UK3hap300_QC_norel.rds")
G <- celiac$genotypes$copy(code = c(0, 1, 2, 0, rep(NA, 252)))
CHR <- celiac$map$chromosome
POS <- celiac$map$physical.pos

system.time(
  corr <- snp_cor(G, ind.col = which(CHR == 22), infos.pos = POS[CHR == 22], ncores = 4)
)
## Before: 12 sec just for chr22 -> 18 sec with just one core
## OpenMP: 6-7 sec with 4 cores and 17 sec with one core
length(corr@x)  # 263,434
object.size(corr) / length(corr@x)  # 12 = 4 for int + 8 for double
length(corr@x) / length(corr)  # ~1%

# Also, the filling of the matrix does not affect performance
# -> Armadillo if now optimized for this and I don't need batch insertions
