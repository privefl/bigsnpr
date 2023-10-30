library(bigsnpr)
NCORES <- nb_cores()

bedfile <- download_1000G("../datasets")
fam2 <- bigreadr::fread2(sub_bed(bedfile, ".fam2"))
ind_eur <- which(fam2$`Super Population` == "EUR")

obj.bed <- bed(bedfile)
obj.svd0 <- bed_autoSVD(obj.bed, ind.row = ind_eur, ncores = NCORES)


#### alternative approach ####

# First clumping
ind_keep <- bed_clumping(obj.bed, ind.row = ind_eur,
                         thr.r2 = 0.5, size = 500, ncores = NCORES)



min.mac <- 10
if (min.mac > 0) {
  mac <- bed_MAF(obj.bed, ind_eur, ind_keep, ncores = NCORES)$mac
  mac.nok <- (mac < min.mac)
  ind_keep <- ind_keep[!mac.nok]
}


obj.svd <- bed_randomSVD(obj.bed,
                         fun.scaling = bed_scaleBinom,
                         ind.row = ind_eur,
                         ind.col = ind_keep,
                         k = 10,
                         ncores = NCORES)

plot(obj.svd, type = "loadings", loadings = 1:10, coeff = 0.4)

obj.pcadapt <- bed_pcadapt(obj.bed, obj.svd$u, ind.row = ind_eur,
                           ind.col = ind_keep, ncores = NCORES)
lpval <- -predict(obj.pcadapt)
pval2 <- 10^attr(obj.pcadapt, "predict")(obj.pcadapt$score)
hist(pval2)

ind_top <- ind_keep[pval2 < 5e-8]
S <- rep(NA, ncol(obj.bed)); S[ind_keep] <- -pval2
ind_top_indep <- bed_clumping(obj.bed, ind.row = ind_eur, S = S,
                              exclude = setdiff(cols_along(obj.bed), ind_top),
                              thr.r2 = 0.05, size = Inf, ncores = NCORES)
(ind_rm <- setdiff(ind_top, ind_top_indep))
stopifnot(all(ind_rm %in% ind_keep))
ind_keep <- setdiff(ind_keep, ind_rm)

options(max.print = 200)
round(100 * cor(obj.svd0$u, obj.svd$u), 1)
