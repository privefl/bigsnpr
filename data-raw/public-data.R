# https://www.dropbox.com/s/k9ptc4kep9hmvz5/1kg_phase1_all.tar.gz?dl=1

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
ped <- data.table::fread("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped", data.table = FALSE)
fam <- dplyr::left_join(fam, ped, by = c("sample.ID" = "Individual ID"))
pop <- data.table::fread("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20131219.populations.tsv", data.table = FALSE)
fam <- dplyr::left_join(fam, pop, by = c("Population" = "Population Code"))

snp$fam$family.ID <- paste(fam$`Super Population`, fam$Population, sep = "_")
snp$fam$paternal.ID <- snp$fam$maternal.ID <- 0L

ind <- which(fam$Relationship == "unrel")
maf <- snp_MAF(G, ind)

bed <- snp_writeBed(snp, tempfile(fileext = ".bed"),
                    ind.row = ind, ind.col = which(maf > 0.05))

rds <- snp_readBed(bed)
snp <- snp_attach(rds)
G <- snp$genotypes
set.seed(1)
pheno <- pkg.paper.PRS::get_pheno(G, 0.8, 10)
snp$fam$affection <- pheno$pheno

snp_writeBed(snp, "tmp-data/public-data.bed")
saveRDS(pheno, file = "tmp-data/public-data-pheno.rds")
