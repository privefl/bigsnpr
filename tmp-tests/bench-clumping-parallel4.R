library(bigsnpr)

obj.bed <- bed("../Dubois2010_data/FinnuncorrNLITUK1UK3hap300_QC_norel.bed")

system.time(
  keep <- bed_clumping(obj.bed, ncores = 4)
)
# Before: 1 -> 62 sec // 4 -> 34 sec
# After:  1 -> 67 sec // 4 -> 50-54 sec
# 5 seconds are always taken by colstats
length(keep)
# Before: 94,777
# After:  94,777
