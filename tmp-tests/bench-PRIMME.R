library(bigsnpr)
obj.bed <- bed("../POPRES_data/POPRES_allchr_QC_norel.bed")
ind.keep <- bed_clumping(obj.bed, thr.r2 = 0.2, size = 1000,
                         exclude = which(obj.bed$map$chromosome %in% c(2:6)),
                         ncores = nb_cores())


system.time(obj.svd <- bed_randomSVD(obj.bed, ind.col = ind.keep))
# 5 sec -> 28 sec

system.time(obj.svd <- bed_autoSVD(
  obj.bed, ind.col = ind.keep, thr.r2 = NA, tukey.coef = 1.5))
# 17.5 sec with no prior       -> 105 sec
# 20 sec while providing both  -> 128 sec
# 22 sec while providing u0    -> 121 sec
# 22 sec while providing v0    -> 111 sec
