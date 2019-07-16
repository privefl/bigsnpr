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
rsid.hapmap3 <- bigreadr::fread2(input = "ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/hapmap3/plink_format/draft_2/hapmap3_r2_b36_fwd.consensus.qc.poly.map.bz2",
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

download_plink(".")
snp_plinkIBDQC("./plink",
               "1000G_phase3_common_hapmap.bed",
               pi.hat = 0.25,
               pruning.args = c(2000, 0.1),
               ncores = NCORES,
               extra.options = "--memory 8000")


zip("1000G_phase3_common_hapmap.zip",
    paste0("1000G_phase3_common_hapmap", c(".bed", ".bim", ".fam")))


system("./plink2 --bfile 1000G_phase3_common_hapmap --freq --out 1000G_phase3_common_hapmap")
af <- bigreadr::fread2("1000G_phase3_common_hapmap.afreq", select = "ALT_FREQS")[[1]]
summary(af)
sum(af > 0.05)

setwd(wd)
