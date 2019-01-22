# https://www.dropbox.com/s/k9ptc4kep9hmvz5/1kg_phase1_all.tar.gz?dl=1

library(bigsnpr)
plink <- download_plink("tmp-data/")
bed <- snp_plinkQC(plink, prefix.in = "tmp-data/1kg_phase1_all",
                   geno = 0, maf = 0.05, hwe = 1e-10, autosome.only = TRUE,
                   extra.options = " --thin 0.05")

bed2 <- snp_plinkIBDQC(plink, bed, pi.hat = 0.05, ncores = nb_cores())

rds <- snp_readBed(bed2)
snp <- snp_attach(rds)
G <- snp$genotypes
counts <- big_counts(G)
sum(counts[4, ])
counts[, 1:10]

fam <- snp$fam
ped <- data.table::fread("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped", data.table = FALSE)
fam <- dplyr::left_join(fam, ped, by = c("sample.ID" = "Individual ID"))
pop <- data.table::fread("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20131219.populations.tsv", data.table = FALSE)
fam <- dplyr::left_join(fam, pop, by = c("Population" = "Population Code"))

snp$fam$family.ID <- paste(fam$`Super Population`, fam$Population, sep = "_")
snp$fam$paternal.ID <- snp$fam$maternal.ID <- 0L

maf <- snp_MAF(G)

snp_writeBed(snp, "tmp-data/public-data.bed", ind.col = which(maf > 0.05))
