G2 <- big_copy(G)
G2[sample(length(G2), 0.1 * length(G2))] <- 3  # NAs
G2[, 1]
snp2 <- snp
snp2$genotypes <- G2
bedfile2 <- snp_writeBed(snp2, tempfile(fileext = ".bed"))

coef <- sapply(cols_along(G2), function(j) {
  summary(lm(U[, 1] ~ G2[, j]))$coefficients[2, 3]
})

obj.bed <- bed(bedfile2)
stats <- bigsnpr:::bed_stats(obj.bed, ind.row, ind.col)
center <- stats$sum / stats$nb_nona_col
scale <- rep(1, length(center))
test <- bigsnpr:::multLinReg(obj.bed, rows_along(G), cols_along(G),
                             center, scale, U[, 1, drop = FALSE])
all.equal(test[, 1], coef)
plot(test[, 1], coef)
mean(test[, 1] / coef)

bed4pcadapt <- pcadapt::read.pcadapt(bedfile2, type = "bed")
obj.pcadapt <- pcadapt::pcadapt(bed4pcadapt, K = 10, min.maf = 0)
test2 <- obj.pcadapt$zscores[, 1]
plot(test2, -coef)
all.equal(test2, -coef)

plink <- download_plink("tmp-data")
prefix_bed <- sub_bed(bedfile2)
# GWAS (linear)
tmp <- tempfile(fileext = ".phe")
bigsnpr:::write.table2(cbind(snp$fam[, 1:2], U), tmp)
system.time(
  system(glue::glue("{plink} --bfile {prefix_bed} --out {prefix_bed}",
                    " --linear --pheno {tmp}",
                    # " --mpheno 1",
                    " --allow-no-sex",
                    " --threads {NCORES}"))
)
gwas2.lin <- bigreadr::fread2(paste0(prefix_bed, ".assoc.linear"))
stat <- ifelse(gwas2.lin$A1 == snp2$map$allele1, gwas2.lin$STAT, -gwas2.lin$STAT)
plot(stat, coef)
plot(stat, -test2)
all.equal(stat, coef)
all.equal(stat, -test2)
