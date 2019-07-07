# source('tmp-tests/test-pcadapt4.R', echo = TRUE)
G2 <- big_copy(G); M <- ncol(G); N <- nrow(G)
nbNA <- VGAM::rbetabinom.ab(M, size = N, shape1 = 0.6, shape2 = 2)
sum(nbNA) / length(G)
indNA <- cbind(
  unlist(
    lapply(nbNA, function(nb) {
      `if`(nb > 0, sample(N, size = nb), NULL)
    })
  ),
  rep(cols_along(G), nbNA)
)
G2[indNA] <- 3
big_counts(G2)[4, ] / nrow(G2)
# G2[sample(length(G2), 0.4 * length(G2))] <- 3  # NAs
G2[, 1]
snp2 <- snp
snp2$genotypes <- G2
bedfile2 <- snp_writeBed(snp2, tempfile(fileext = ".bed"))

ind.keep2 <- bed_clumping(bedfile2)
svd2 <- bed_randomSVD(bed(bedfile2), ind.col = ind.keep2)
U2 <- svd2$u
plot(U, U2)
plot(U[, 1:2], U2[, 1:2], pch = 20); abline(0, 1, col = "red")

coef <- sapply(cols_along(G2), function(j) {
  summary(lm(U2[, 9] ~ G2[, j]))$coefficients[2, 3]
})
# all.equal(coef, big_univLinReg(G, U[, 9])$score)

obj.bed <- bed(bedfile2)
test <- bigsnpr:::multLinReg(obj.bed, rows_along(G), cols_along(G), U2)
all.equal(test[, 9], coef)

bed4pcadapt <- pcadapt::read.pcadapt(bedfile2, type = "bed")
obj.pcadapt <- pcadapt::pcadapt(bed4pcadapt, K = 10, min.maf = 0)
test2 <- -obj.pcadapt$zscores[, 1]
plot(test2, test[, 1])
plot(test2, big_univLinReg(G, U[, 1])$score)
all.equal(test2, big_univLinReg(G, U[, 1])$score)

# plink <- download_plink("tmp-data")
# prefix_bed <- sub_bed(bedfile2)
# # GWAS (linear)
# tmp <- tempfile(fileext = ".phe")
# bigsnpr:::write.table2(cbind(snp$fam[, 1:2], U2), tmp)
# system.time(
#   system(glue::glue("{plink} --bfile {prefix_bed} --out {prefix_bed}",
#                     " --linear --pheno {tmp}",
#                     # " --mpheno 1",
#                     " --allow-no-sex",
#                     " --threads 2"))
# )
# gwas2.lin <- bigreadr::fread2(paste0(prefix_bed, ".assoc.linear"))
# stat <- ifelse(gwas2.lin$A1 == snp2$map$allele1, gwas2.lin$STAT, -gwas2.lin$STAT)
# plot(stat, coef)
# plot(stat, test2)
# all.equal(stat, coef)
# all.equal(stat, test2)
c(
  # PLINK = all.equal(stat, big_univLinReg(G, U[, 1])$score),
  bigsnpr = all.equal(test[, 1], big_univLinReg(G, U[, 1])$score),
  pcadapt = all.equal(test2, big_univLinReg(G, U[, 1])$score)
)
k <- 1
c(
  bigsnpr = all.equal(test[, k], big_univLinReg(G, U[, k])$score),
  pcadapt = all.equal(-obj.pcadapt$zscores[, k], big_univLinReg(G, U[, k])$score)
)


library(ggplot2)
cowplot::plot_grid(
  qplot(test[, 1], big_univLinReg(G, U[, 1])$score, color = big_counts(G2)[4, ]) +
    geom_abline(color = "red") +
    theme_bigstatsr() +
    scale_color_viridis_c() +
    theme(legend.position = c(0.8, 0.3)) +
    labs(color = "#NAs", x = "Stat bigsnpr"),
  qplot(test2, big_univLinReg(G, U[, 1])$score, color = big_counts(G2)[4, ]) +
    geom_abline(color = "red") +
    theme_bigstatsr() +
    scale_color_viridis_c() +
    theme(legend.position = c(0.8, 0.3)) +
    labs(color = "#NAs", x = "Stat pcadapt")
)
