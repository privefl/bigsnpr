setwd("tmp-data/")

# https://cran.r-project.org/web/packages/plinkQC/vignettes/Genomes1000.pdf
download.file("https://www.dropbox.com/s/afvvf1e15gqzsqo/all_phase3.pgen.zst?dl=1",
              destfile = "all_phase3.pgen.zst")
download.file("https://www.dropbox.com/s/op9osq6luy3pjg8/all_phase3.pvar.zst?dl=1",
              destfile = "all_phase3.pvar.zst")
download.file("https://www.dropbox.com/s/yozrzsdrwqej63q/phase3_corrected.psam?dl=1",
              destfile = "all_phase3.psam")

library(bigsnpr)
plink2 <- download_plink2(".")

system("./plink2 --zst-decompress all_phase3.pgen.zst > all_phase3.pgen")
system("./plink2 --zst-decompress all_phase3.pvar.zst > all_phase3.pvar")
system("./plink2 --pfile all_phase3 --make-bed --max-alleles 2 --out all_phase3")

bed <- snp_plinkQC("./plink2", "all_phase3", geno = 0, hwe = 1e-10,
                   autosome.only = TRUE, extra.options = "--thin 0.05")
# snp_plinkIBDQC(download_plink(), bed, ncores = nb_cores())

snp_readBed("all_phase3_QC.bed")
snp <- snp_attach("all_phase3_QC.rds")
G <- snp$genotypes
counts <- big_counts(G)
sum(counts[4, ])
counts[, 1:10]

pheno <- pkg.paper.PRS::get_pheno(G, 0.5, 1000)
snp$fam$affection <- pheno$pheno

snp_writeBed(snp, "public-data2.bed")
saveRDS(pheno, "public-data2-pheno.rds")

system.time(
  gwas <- big_univLogReg(G, pheno$pheno, ncores = nb_cores())
)
plot(gwas)
snp_qq(gwas)
