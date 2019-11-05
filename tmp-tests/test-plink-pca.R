library(bigsnpr)

bedfile <- download_1000G("tmp-data/")
obj.bed <- bed(bedfile)
CHR <- obj.bed$map$chromosome
POS <- obj.bed$map$physical.pos
ind.keep <- bed_clumping(obj.bed, ncores = nb_cores(),
                         exclude = snp_indLRLDR(CHR, POS))

tmp <- tempfile()
write(obj.bed$map$marker.ID[ind.keep], tmp, ncolumns = 1)
plink2 <- download_plink2("tmp-data/")
system.time(
  system(glue::glue(
    "{plink2} --pca approx",
    " --bfile {sub_bed(bedfile)}",
    " --extract {tmp}",
    " --out {tmp}",
    " --memory 6000"
  ))
) # 66 sec -> 490 -> 90

val <- bigreadr::fread2(paste0(tmp, ".eigenval"))
vec <- bigreadr::fread2(paste0(tmp, ".eigenvec"))
plot(vec[3:4 + 6])

system.time(
  svd <- bed_randomSVD(obj.bed, ind.col = ind.keep, ncores = nb_cores())
) # 37 sec -> 320 -> 99


plot(as.matrix(vec[-(1:2)]), svd$u, pch = 20)
abline(0, 1, col = "red"); abline(0, -1, col = "red")
round(100 * cor(as.matrix(vec[-(1:2)]), svd$u), 1)
