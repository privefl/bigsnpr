# download.file("https://www.dropbox.com/s/k9ptc4kep9hmvz5/1kg_phase1_all.tar.gz?dl=1",
#               destfile = "tmp-data/1kg_phase1_all.tar.gz")
# untar("tmp-data/1kg_phase1_all.tar.gz", exdir = "tmp-data/")

library(bigsnpr)
plink <- download_plink("tmp-data/")
bed <- snp_plinkQC(plink, prefix.in = "tmp-data/1kg_phase1_all",
                   geno = 0, maf = 0.05, hwe = 1e-10,
                   extra.options = " --chr 2,6,8 --thin 0.1")

rds <- snp_readBed(bed)
snp <- snp_attach("tmp-data/1kg_phase1_all_QC.rds")
G <- snp$genotypes
counts <- big_counts(G)
sum(counts[4, ])
counts[, 1:10]

fam <- snp$fam
ped <- bigreadr::fread2("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped")
fam <- dplyr::left_join(fam, ped, by = c("sample.ID" = "Individual ID"))
pop <- bigreadr::fread2("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20131219.populations.tsv")
fam <- dplyr::left_join(fam, pop, by = c("Population" = "Population Code"))

snp$fam$family.ID <- paste(fam$`Super Population`, fam$Population, sep = "_")
snp$fam$paternal.ID <- snp$fam$maternal.ID <- 0L

ind_norel <- which(fam$Relationship == "unrel")
maf <- snp_MAF(G, ind_norel)

bed <- snp_writeBed(snp, tempfile(fileext = ".bed"),
                    ind.row = ind_norel, ind.col = which(maf > 0.05))

rds <- snp_readBed(bed)
snp <- snp_attach(rds)
G <- snp$genotypes
CHR <- snp$map$chromosome
POS <- snp$map$physical.pos
set.seed(1)
# devtools::install_github("privefl/paper2-PRS/pkg.paper.PRS")
pheno <- pkg.paper.PRS::get_pheno(G, 0.8, 100)
snp$fam$affection <- pheno$pheno + 1
s <- scale(G[, pheno$set]) %*% pheno$effects + rnorm(nrow(G), sd = 2)
obj.svd <- snp_autoSVD(G, CHR, POS, thr.r2 = 0.1)
plot(obj.svd, type = "scores")
gwas <- big_univLinReg(G, s, covar.train = obj.svd$u, ncores = nb_cores())
plot(gwas)
snp_manhattan(gwas, CHR, POS)
y_hat <- mean(pheno$pheno)
plot(pheno$effects, gwas$estim[pheno$set] * y_hat * (1 - y_hat), pch = 20); abline(0, 1, col = "red")
sumstats <- cbind.data.frame(snp$map[-3], beta = gwas$estim, p = predict(gwas, log10 = FALSE))

snp_writeBed(snp, "tmp-data/public-data.bed")
saveRDS(pheno, file = "tmp-data/public-data-pheno.rds")
bigreadr::fwrite2(sumstats, "tmp-data/public-data-sumstats.txt")
zip("data-raw/public-data.zip",
    paste0("tmp-data/public-data", c(".bed", ".bim", ".fam", "-pheno.rds", "-sumstats.txt")))
