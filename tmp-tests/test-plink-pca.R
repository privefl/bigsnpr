library(bigsnpr)

obj.bed <- bed("tmp-data/1000G_phase3_common_hapmap_norel.bed")
plink2 <- download_plink2("tmp-data/")
system.time(
  system(glue::glue(
    "{plink2} --pca approx",
    " --bfile tmp-data/1000G_phase3_common_hapmap_norel",
    " --out tmp-data/test-plink"
  ))
) # 66 sec

val <- bigreadr::fread2("tmp-data/test-plink.eigenval")
vec <- bigreadr::fread2("tmp-data/test-plink.eigenvec")
plot(vec[3:4 + 6])

system.time(
  svd <- bed_randomSVD(obj.bed, ncores = nb_cores())
) # 37 sec


plot(as.matrix(vec[-(1:2)]), svd$u, pch = 20)
abline(0, 1, col = "red"); abline(0, -1, col = "red")
round(100 * cor(as.matrix(vec[-(1:2)]), svd$u), 1)
