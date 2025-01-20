library(bigsnpr)

celiac <- snp_attach("../Dubois2010_data/FinnuncorrNLITUK1UK3hap300_QC_norel.rds")
dim(G <- celiac$genotypes)  # 15155 x 281122
CHR <- celiac$map$chromosome
POS <- celiac$map$physical.pos
y <- celiac$fam$affection - 1
NCORES <- nb_cores()

pop <- snp_getSampleInfos(celiac, "../Dubois2010_data/FinnNLITUK1UK3.clusterv2")[[1]]
pop2 <- c("Netherlands", "Italy", "UK1", "UK2", "Finland")[pop]
table(pop2, exclude = NULL)
# Finland       Italy Netherlands         UK1         UK2
#    2436        1035        1623        3325        6736

big_counts(G, ind.col = 1:10)
G2 <- snp_fastImputeSimple(G, method = "mean2", ncores = NCORES)

obj.svd <- snp_autoSVD(G2, CHR, POS, ncores = NCORES)
plot(obj.svd)
plot(obj.svd, type = "scores", scores = 1:6, coef = 0.6)
PC <- predict(obj.svd)[, 1:6]

pop_centers <- bigutilsr::geometric_median(PC, by_grp = pop2)

ldist_to_centers <- apply(pop_centers, 1, function(center) {
  log(rowSums(sweep(PC, 2, center, '-')^2))
})
hist(ldist_to_centers[pop2 != "Finland", "Finland"])
hist(ldist_to_centers[pop2 == "Italy", "Italy"])
hist(ldist_to_centers[pop2 != "Italy", "Italy"])
hist(ldist_to_centers[pop2 != "UK2", "UK2"])

ind_NW <- which(ldist_to_centers[, "UK2"] < 6)
ind_Fin <- which(ldist_to_centers[, "Finland"] < 7)
ind_It <- which(ldist_to_centers[, "Italy"] < 6)

run_GWAS <- function(ind) {
  big_univLinReg(G2, y.train = y[ind], ind.train = ind,
                 covar.train = cbind(celiac$fam$sex, PC)[ind, ],
                 ncores = NCORES)
}

gwas_NW <- run_GWAS(ind_NW)
gwas_Fin <- run_GWAS(ind_Fin)
gwas_It <- run_GWAS(ind_It)

plot_grid(
  snp_manhattan(gwas_NW,  CHR, POS, npoints = 20e3) +
    ggplot2::scale_y_log10(),
  snp_manhattan(gwas_Fin, CHR, POS, npoints = 20e3) +
    ggplot2::scale_y_log10(),
  snp_manhattan(gwas_It,  CHR, POS, npoints = 20e3) +
    ggplot2::scale_y_log10(),
  ncol = 1
)

ind_sig <- which(predict(gwas_NW, log10 = FALSE) < 1e-20)
plot(gwas_NW$estim[ind_sig], gwas_Fin$estim[ind_sig]); abline(0, 1, col = "red")
plot(gwas_NW$estim[ind_sig], gwas_It$estim[ind_sig]); abline(0, 1, col = "red")

reg <- deming::deming(gwas_NW$estim[ind_sig] ~ gwas_Fin$estim[ind_sig] + 0,
                      xstd = gwas_Fin$std.err[ind_sig], ystd = gwas_NW$std.err[ind_sig])
reg

deming::deming(gwas_NW$estim[ind_sig] ~ gwas_It$estim[ind_sig] + 0,
               xstd = gwas_It$std.err[ind_sig], ystd = gwas_NW$std.err[ind_sig])


ind_chr6 <- which(CHR == 6)
POS2 <- snp_asGeneticPos(CHR[ind_chr6], POS[ind_chr6])
corr_chr6_NW <- snp_cor(G, ind.row = ind_NW, ind.col = ind_chr6,
                        infos.pos = POS2, size = 3 / 1000, ncores = NCORES)
corr_chr6_Fin <- snp_cor(G, ind.row = ind_Fin, ind.col = ind_chr6,
                         infos.pos = POS2, size = 3 / 1000, ncores = NCORES)
corr_chr6_It <- snp_cor(G, ind.row = ind_It, ind.col = ind_chr6,
                        infos.pos = POS2, size = 3 / 1000, ncores = NCORES)

ind <- Matrix::which(corr_chr6_NW^2 > 0.1 & corr_chr6_Fin^2 > 0.1, arr.ind = TRUE)
ind2 <- ind[ind[, 1] > ind[, 2], ]
plot(corr_chr6_Fin[ind2], corr_chr6_NW[ind2], pch = 20); abline(0, 1, col = "red")


ind3 <- Matrix::which(corr_chr6_NW^2 > 0.1 & corr_chr6_It^2 > 0.1, arr.ind = TRUE)
ind4 <- ind3[ind3[, 1] > ind3[, 2], ]
plot(corr_chr6_It[ind4], corr_chr6_NW[ind4], pch = 20); abline(0, 1, col = "red")

plot(corr_chr6_It[ind4], corr_chr6_Fin[ind4], pch = 20); abline(0, 1, col = "red")


ind_UK2 <- which(ldist_to_centers[, "UK2"] < 5 & pop2 == "UK2")

gwas_UK1 <- run_GWAS(ind_UK1)
gwas_UK2 <- run_GWAS(ind_UK2)

plot_grid(
  snp_manhattan(gwas_UK1,  CHR, POS, npoints = 20e3) +
    ggplot2::scale_y_log10(),
  snp_manhattan(gwas_UK2, CHR, POS, npoints = 20e3) +
    ggplot2::scale_y_log10(),
  ncol = 1
)
corr_chr6_UK2 <- snp_cor(G, ind.row = ind_UK2, ind.col = ind_chr6,
                         infos.pos = POS2, size = 3 / 1000, ncores = NCORES)
plot(corr_chr6_UK1[ind2], corr_chr6_UK2[ind2], pch = 20); abline(0, 1, col = "red")


ind_Net <- which(ldist_to_centers[, "Netherlands"] < 5 & pop2 == "Netherlands")
corr_chr6_Net <- snp_cor(G, ind.row = ind_Net, ind.col = ind_chr6,
                         infos.pos = POS2, size = 3 / 1000, ncores = NCORES)
plot(corr_chr6_Net[ind2], corr_chr6_UK2[ind2], pch = 20); abline(0, 1, col = "red")


ind_UK1 <- sample(which(ldist_to_centers[, "UK1"] < 5 & pop2 == "UK1"), length(ind_Net))
corr_chr6_UK1 <- snp_cor(G, ind.row = ind_UK1, ind.col = ind_chr6,
                         infos.pos = POS2, size = 3 / 1000, ncores = NCORES)
plot(corr_chr6_UK1[ind2], corr_chr6_UK2[ind2], pch = 20); abline(0, 1, col = "red")


ind_max <- which.min(predict(gwas_NW))
ind5 <- ind_max + c(-4:2, 4, 10)
plot_grid(
  snp_manhattan(gwas_NW[ind5, ],  CHR[ind5], POS[ind5]) +
    ggplot2::scale_y_log10(),
  snp_manhattan(gwas_Fin[ind5, ], CHR[ind5], POS[ind5]) +
    ggplot2::scale_y_log10(),
  snp_manhattan(gwas_It[ind5, ],  CHR[ind5], POS[ind5]) +
    ggplot2::scale_y_log10(),
  ncol = 1
)

ind6 <- match(ind5, ind_chr6)
options(max.print = 200)
rbind(
  corr_chr6_NW[ind6, ind6[5]],
  corr_chr6_Fin[ind6, ind6[5]],
  corr_chr6_It[ind6, ind6[5]])
