library(bigsnpr)
NCORES <- nb_cores()

bedfile <- "../datasets/pcadaptissue83/pcadapt_input_for_prive/input.plink.bed"

obj.bed <- bed(bedfile)
# obj.svd0 <- bed_autoSVD(obj.bed, ncores = NCORES)


#### alternative approach ####

# First clumping
# ind_keep <- bed_clumping(obj.bed, thr.r2 = 0.2, size = 500, ncores = NCORES)

obj.bed$.map$chromosome <- 1

min.mac <- 10
mac <- bed_MAF(obj.bed, ncores = NCORES)$mac
ind_keep <- which(mac > min.mac)

obj.svd0 <- bed_randomSVD(obj.bed,
                         fun.scaling = bed_scaleBinom,
                         k = 20,
                         ncores = NCORES)
plot(obj.svd0)
plot(obj.svd0, type = "scores", scores = 9:10 + 8) +
  ggplot2::aes(color = as.factor(pop$V1))

obj.pcadapt0_full <- bed_pcadapt(obj.bed, obj.svd0$u, ncores = NCORES)
plot(obj.pcadapt0_full)
plot(obj.pcadapt0_full, "Manhattan") +
  ggplot2::geom_hline(yintercept = -log10(5e-8), color = "red")
sum(predict(obj.pcadapt0_full, log10 = FALSE) < 5e-8)  # 459


#### REPEAT UNTIL NO IND_RM ####

obj.svd <- bed_randomSVD(obj.bed,
                         fun.scaling = bed_scaleBinom,
                         ind.col = ind_keep,
                         k = 14,
                         ncores = NCORES)

plot(obj.svd)
# plot(obj.svd, type = "loadings", loadings = 1:10, coeff = 0.4)

obj.pcadapt <- bed_pcadapt(obj.bed, obj.svd$u, ind.col = ind_keep, ncores = NCORES)
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

#### NO IND_RM ####


plot(obj.svd, type = "scores", scores = 1:8, coeff = 0.5)



obj.pcadapt_full <- bed_pcadapt(obj.bed, obj.svd$u, ncores = NCORES)
plot(obj.pcadapt_full)
plot(obj.pcadapt_full, "Manhattan") +
  ggplot2::geom_hline(yintercept = -log10(5e-8), color = "red")
sum(predict(obj.pcadapt_full, log10 = FALSE) < 5e-8)  # 394
pop <- bigreadr::fread2("../datasets/pcadaptissue83/pcadapt_input_for_prive/pop_ind.txt",
                        header = FALSE)
plot(obj.svd, type = "scores", scores = 9:10 + 4) +
  ggplot2::aes(color = as.factor(pop$V1))
