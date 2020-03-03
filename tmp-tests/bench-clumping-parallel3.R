library(bigsnpr)

celiac <- snp_attach("../Dubois2010_data/FinnuncorrNLITUK1UK3hap300_QC_norel.rds")
G <- celiac$genotypes$copy(code = c(0, 1, 2, 0, rep(NA, 252)))
CHR <- celiac$map$chromosome
POS <- celiac$map$physical.pos

ind.row = rep(rows_along(G), 5)
ind.col <- which(CHR %in% 1:4)
G2 <- big_copy(G, ind.row, ind.col)

RcppParallel::setThreadOptions(2)
system.time(
  keep <- snp_clumping(G2, infos.chr = CHR[ind.col], infos.pos = POS[ind.col],
                       ncores = 2, thr.r2 = 0.95)
)
# Before: 1 -> 55 sec // 4 -> 22 sec // r2=0.95 -> 31
# After:  1 -> 58 sec // 4 -> 35 sec // r2=0.95 -> 40 vs 78 // r2=0.01 -> 38 vs 66
# 4 seconds are always taken by colstats
length(keep)
# Before: 27,067
# After:  27,067
