library(bigsnpr)
NCORES <- nb_cores()

bedfile <- "../Dubois2010_data/FinnuncorrNLITUK1UK3hap300_QC_norel.bed"

obj.bed <- bed(bedfile)
obj.svd0 <- bed_autoSVD(obj.bed, ncores = NCORES)


#### alternative approach ####

# First clumping
ind_keep <- bed_clumping(obj.bed, thr.r2 = 0.2, size = 500, ncores = NCORES)



min.mac <- 10
if (min.mac > 0) {
  mac <- bed_MAF(obj.bed, ind.col = ind_keep, ncores = NCORES)$mac
  C <- (mac < min.mac)
  ind_keep <- ind_keep[!mac.nok]
}


obj.svd <- bed_randomSVD(obj.bed,
                         fun.scaling = bed_scaleBinom,
                         ind.col = ind_keep,
                         k = 10,
                         ncores = NCORES)

plot(obj.svd, type = "loadings", loadings = 1:10, coeff = 0.4)

obj.pcadapt <- bed_pcadapt(obj.bed, obj.svd$u,
                           ind.col = ind_keep, ncores = NCORES)
pval2 <- 10^attr(obj.pcadapt, "predict")(obj.pcadapt$score)
hist(pval2)

ind_top <- ind_keep[pval2 < 5e-8]
S <- rep(NA, ncol(obj.bed)); S[ind_keep] <- -pval2
ind_top_indep <- bed_clumping(obj.bed, S = S,
                              exclude = setdiff(cols_along(obj.bed), ind_top),
                              thr.r2 = 0.05, size = Inf, ncores = NCORES)
(ind_rm <- setdiff(ind_top, ind_top_indep))
stopifnot(all(ind_rm %in% ind_keep))
ind_keep <- setdiff(ind_keep, ind_rm)

options(max.print = 200)
round(100 * cor(obj.svd0$u, obj.svd$u), 1)

plot(obj.svd0, type = "scores", scores = 1:8, coeff = 0.5)
