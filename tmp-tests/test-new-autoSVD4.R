library(bigsnpr)
NCORES <- nb_cores()

bedfile <- "../datasets/pcadaptissue83/pcadapt_input_for_prive/input.plink.bed"
mat <- pcadapt::bed2matrix(bedfile)
# Verif
hist(colMeans(mat, na.rm = TRUE))
hist(colMeans(is.na(mat)))
hist(rowMeans(is.na(mat)))


pop <- bigreadr::fread2("../datasets/pcadaptissue83/pcadapt_input_for_prive/pop_ind.txt",
                        header = FALSE)

# adding some null SNPs
N <- nrow(mat)
M <- 50e3
mat_null <- sapply(runif(M, min = 0.05, max = 0.5), function(af) {
  rbinom(N, size = 2, prob = af)
})
mat2 <- cbind(mat, mat_null)
pcadapt::read.pcadapt(mat2, type = "lfmm")

fake <- snp_fake(N, ncol(mat2))
mat2[is.na(mat2)] <- 3L
fake$genotypes[] <- mat2
tmp <- tempfile(fileext = ".bed")
snp_writeBed(fake, tmp)

obj.bed <- bed(tmp)
# obj.svd0 <- bed_autoSVD(obj.bed, ncores = NCORES)


#### alternative approach ####


min.mac <- 10
mac <- bed_MAF(obj.bed, ncores = NCORES)$mac
ind_keep <- which(mac > min.mac)

obj.svd0 <- bed_randomSVD(obj.bed,
                          fun.scaling = bed_scaleBinom,
                          k = 20,
                          ncores = NCORES)
plot(obj.svd0)
plot(obj.svd0, type = "scores", scores = 9:20) +
  ggplot2::aes(color = as.factor(pop$V1))
plot(obj.svd0, type = "scores", scores = 11:12) +
  ggplot2::aes(color = as.factor(pop$V1))

obj.pcadapt0_full <- bed_pcadapt(obj.bed, obj.svd0$u[, 1:12], ncores = NCORES)
plot(obj.pcadapt0_full)
plot(obj.pcadapt0_full, "Manhattan") +
  ggplot2::geom_hline(yintercept = -log10(5e-8), color = "red")
sum(predict(obj.pcadapt0_full, log10 = FALSE) < 5e-8)  # 459


#### REPEAT UNTIL NO IND_RM ####

obj.svd <- bed_randomSVD(obj.bed,
                         fun.scaling = bed_scaleBinom,
                         ind.col = ind_keep,
                         k = 12,
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
plot(obj.svd, type = "scores", scores = 5:6, coeff = 0.5) +
  ggplot2::aes(color = as.factor(pop$V1))
plot(obj.svd, type = "scores", scores = 7:8, coeff = 0.5) +
  ggplot2::aes(color = as.factor(pop$V1))
plot(obj.svd, type = "scores", scores = 9:10, coeff = 0.5) +
  ggplot2::aes(color = as.factor(pop$V1))
plot(obj.svd, type = "scores", scores = 11:12, coeff = 0.5) +
  ggplot2::aes(color = as.factor(pop$V1))



obj.pcadapt_full <- bed_pcadapt(obj.bed, obj.svd$u[, 1:10], ncores = NCORES)
plot(obj.pcadapt_full)
plot(obj.pcadapt_full, "Manhattan") +
  ggplot2::geom_hline(yintercept = -log10(5e-8), color = "red")
sum(predict(obj.pcadapt_full, log10 = FALSE) < 5e-8)  # 394

ind <- cols_along(mat)
plot(obj.pcadapt0_full$score[ind], obj.pcadapt_full$score[ind])
abline(0, 1, col = "red")
