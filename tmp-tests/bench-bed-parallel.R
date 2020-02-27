#### Search for 'sendData.SOCKnode' in profiler report ####
library(bigsnpr)
bigsnp <- snp_attachExtdata()
N <- ncol(bigsnp$genotypes)
bedfile <- snp_writeBed(bigsnp, bedfile = tempfile(fileext = ".bed"),
                        ind.row = 1:100, ind.col = rep(1:N, each = 2000))

file.size(bedfile) / 1024**2  # 217 MB

obj.bed <- bed(bedfile)
dim(obj.bed)  # 100 x 9,084,000
object.size(obj.bed)  # 672 B
object.size(obj.bed$map) / 1024**2 # 312 MB

dup <- duplicated(obj.bed$map$physical.pos)
ind.excl <- which(dup)
ind.keep <- which(!dup)


system.time(test  <- bed_clumping(obj.bed, exclude = ind.excl))             # 2 sec
system.time(test2 <- bed_clumping(obj.bed, exclude = ind.excl, ncores = 2)) # 14 sec
# 2 / 14 if storing fam/map but not nrow/ncol
# 3-4 / 8-9  if storing nrow/ncol but not fam/map
# 2 / 14 if storing fam/map and nrow/ncol
# 2 / 6 when using bed_light

system.time(test  <- bed_autoSVD(obj.bed, ind.col = ind.keep)) # 3 sec
system.time(test2 <- bed_autoSVD(obj.bed, ind.col = ind.keep, ncores = 2))
# 35 sec -> 21 sec -> 12 sec (but still have to register clusters many times)

system.time(test  <- bed_counts(obj.bed, ind.col = ind.keep)) # 0 sec
system.time(test2 <- bed_counts(obj.bed, ind.col = ind.keep, ncores = 2))
# 10 sec -> 2 sec OK

system.time(test  <- bed_MAF(obj.bed, ind.col = ind.keep)) # 0 sec
system.time(test2 <- bed_MAF(obj.bed, ind.col = ind.keep, ncores = 2))
# 10 sec -> 2 sec OK

mac <- bed_MAF(obj.bed, ind.col = ind.keep)$mac
system.time(U <- bed_randomSVD(obj.bed, ind.col = ind.keep[mac > 5], ncores = 2)$u)
system.time(test  <- bed_pcadapt(obj.bed, U.row = U, ind.col = ind.keep[-1469])) # 0 sec
system.time(test2 <- bed_pcadapt(obj.bed, U.row = U, ind.col = ind.keep[-1469], ncores = 2))
# 10 sec -> 2 sec OK
