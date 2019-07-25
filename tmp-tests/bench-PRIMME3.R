library(bigsnpr)
obj.bed <- bed("../POPRES_data/POPRES_allchr_QC_norel.bed")
ind.keep <- bed_clumping(obj.bed, thr.r2 = 0.2, size = 1000,
                         ncores = nb_cores())
system.time(obj.svd <- bed_randomSVD(obj.bed, ind.col = ind.keep, verbose = TRUE))
# Before: 57 sec
# After:  51 sec
plot(obj.svd)
plot(obj.svd, type = "scores")

u0 <- obj.svd$u
v0 <- obj.svd$v
plot(u0, u0 + rnorm(length(u0), sd = 0.01), pch = 20)
# plot(v0, v0 + rnorm(length(v0), sd = 0.01), pch = 20)
system.time(obj.svd2 <- bed_randomSVD(
  obj.bed, ind.col = ind.keep, verbose = TRUE, u0 = u0, v0 = v0)) # 5 sec
system.time(obj.svd2 <- bed_randomSVD(
  obj.bed, ind.col = ind.keep, verbose = TRUE,
  u0 = u0 + rnorm(length(u0), sd = 0.01))) # 53 sec
system.time(obj.svd2 <- bed_randomSVD(
  obj.bed, ind.col = ind.keep, verbose = TRUE,
  u0 = u0 + rnorm(length(u0), sd = 0.01), v0 = v0)) # 58 sec
system.time(obj.svd2 <- bed_randomSVD(
  obj.bed, ind.col = ind.keep, verbose = TRUE,
  v0 = v0 + rnorm(length(v0), sd = 0.01))) # 56 sec
system.time(obj.svd2 <- bed_randomSVD(
  obj.bed, ind.col = ind.keep, verbose = TRUE,
  u0 = u0, v0 = v0 + rnorm(length(v0), sd = 0.01))) # 5 sec
