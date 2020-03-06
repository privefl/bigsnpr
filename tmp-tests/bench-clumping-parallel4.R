library(bigsnpr)

obj.bed <- bed("../Dubois2010_data/FinnuncorrNLITUK1UK3hap300_QC_norel.bed")

system.time(
  keep <- bed_clumping(obj.bed, ncores = 4)
)
# Before: 1 -> 62 sec // 4 -> 34 sec
# After:  1 -> 65 sec // 4 -> 29 sec
length(keep)
# Before: 94,777
# After:  94,777

system.time(print(str(
  bed_prodVec(obj.bed, rep(1, ncol(obj.bed)), ncores = 4)
)))
# 4.8 -> 1.6

system.time(print(str(
  bed_cprodVec(obj.bed, rep(1, nrow(obj.bed)), ncores = 4)
)))
# 4.0 -> 1.4
