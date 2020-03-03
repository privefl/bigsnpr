library(bigsnpr)

celiac <- snp_attach("../Dubois2010_data/FinnuncorrNLITUK1UK3hap300_QC_norel.rds")
G <- celiac$genotypes$copy(code = c(0, 1, 2, 0, rep(NA, 252)))
CHR <- celiac$map$chromosome
POS <- celiac$map$physical.pos

RcppParallel::setThreadOptions(1)
system.time(
  keep <- snp_clumping(G, infos.chr = CHR, infos.pos = POS, ncores = 4)
)
# Before: 1 -> 41 sec // 4 -> 20 sec
# After:  1 -> 43 sec // 4 -> 28 sec
# 3 seconds are always taken by colstats
length(keep)
# Before: 94,831
# After: 94,831

## Previous (2017) number of computations -> 112 sec
# 301064 300555 271914 267613 240766 223591 248161 211148 253087 221697 204203
# 191838 157228 166405 140885 152634 124355 126631 130637 96650 98980 109603
## Current (2019) number of computations -> 63 sec
# 165616 170755 150297 146775 136094 127221 130405 118085 134290 123220 115638
# 106934 89267 89404 80543 87333 71895 75237 75847 60717 52321 60195
## New test sorting by dist (2020) number of computations -> 68 sec
# 166704 171893 151299 147557 136813 127749 131192 118684 135060 123544 116108
# 107556 89840 90108 81171 87401 72247 75456 76036 60677 52476 60435
