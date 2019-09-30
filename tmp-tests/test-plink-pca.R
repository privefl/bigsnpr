library(bigsnpr)

bedfile <- download_1000G("tmp-data/")
plink2 <- download_plink2("tmp-data/")
tmp <- tempfile()
system.time(
  system(glue::glue(
    "{plink2} --pca approx",
    " --bfile {sub_bed(bedfile)}",
    " --out {tmp}",
    " --memory 6000"
  ))
) # 66 sec -> 989.173  43.994 490.083 ??

val <- bigreadr::fread2(paste0(tmp, ".eigenval"))
vec <- bigreadr::fread2(paste0(tmp, ".eigenvec"))
plot(vec[3:4 + 6])

obj.bed <- bed(bedfile)
system.time(
  svd <- bed_randomSVD(obj.bed, ncores = nb_cores())
) # 37 sec -> 320


plot(as.matrix(vec[-(1:2)]), svd$u, pch = 20)
abline(0, 1, col = "red"); abline(0, -1, col = "red")
round(100 * cor(as.matrix(vec[-(1:2)]), svd$u), 1)
