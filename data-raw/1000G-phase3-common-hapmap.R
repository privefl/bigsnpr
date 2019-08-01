wd <- setwd("tmp-data")

# https://cran.r-project.org/web/packages/plinkQC/vignettes/Genomes1000.pdf
download.file("https://www.dropbox.com/s/afvvf1e15gqzsqo/all_phase3.pgen.zst?dl=1",
              destfile = "all_phase3.pgen.zst")
download.file("https://www.dropbox.com/s/0nz9ey756xfocjm/all_phase3.pvar.zst?dl=1",
              destfile = "all_phase3.pvar.zst")
download.file("https://www.dropbox.com/s/yozrzsdrwqej63q/phase3_corrected.psam?dl=1",
              destfile = "all_phase3.psam")

library(bigsnpr)
plink2 <- download_plink2(".")

# Uncompress file
system("./plink2 --zst-decompress all_phase3.pgen.zst > all_phase3.pgen")
system("./plink2 --zst-decompress all_phase3.pvar.zst > all_phase3.pvar")

# Get rsIDs of variants in HapMap3
rsid.hapmap3 <- bigreadr::fread2("ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/hapmap3/plink_format/draft_2/hapmap3_r2_b36_fwd.consensus.qc.poly.map.bz2",
                                 select = 2)
writeLines(rsid.hapmap3[[1]], "hapmap3-rsid.txt")

# Filtering to HapMap3 variants + QC + convertion to bed
bed <- snp_plinkQC(
  plink.path = plink2,
  prefix.in = "all_phase3",
  file.type = "--pfile",
  prefix.out = "1000G_phase3_common_hapmap",
  autosome.only = TRUE,
  extra.options = "--max-alleles 2 --memory 8000 --extract hapmap3-rsid.txt"
)
file.size("1000G_phase3_common_hapmap.bed") / 1024^2  # 787 MB
unlink("all_phase3.*")

# Remove related individuals
bed2 <- snp_plinkKINGQC(plink2, bed, thr.king = 0.177, ncores = nb_cores())

# Get additional information on individuals
fam <- bed(bed2)$fam
ped <- bigreadr::fread2("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped")
fam2 <- dplyr::left_join(fam[c(2, 5)], ped[c(2, 5, 7)],
                         by = c(sample.ID = "Individual ID", sex = "Gender"))
pop <- bigreadr::fread2("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20131219.populations.tsv")
fam3 <- dplyr::left_join(fam2, pop[1:3], by = c("Population" = "Population Code"))
str(fam3)
bigreadr::fwrite2(fam3, sub_bed(bed2, ".fam2"), sep = "\t")
readLines("1000G_phase3_common_hapmap_norel.fam2", n = 5)

# Zip files
zip("1000G_phase3_common_hapmap.zip",
    paste0("1000G_phase3_common_hapmap_norel", c(".bed", ".bim", ".fam", ".fam2")))
file.size("1000G_phase3_common_hapmap.zip") / 1024^2  # 305 MB

unlink(list.files(pattern = "^1000G_phase3_common_hapmap.*[^(\\.zip)]$"))


setwd(wd)
