files <- gsub("\\.ped$", "", list.files("../hapmap3_r2_b36_fwd.qc.poly/",
                                        pattern = "hapmap3_r3_b36_fwd.[A-Z]{3}.qc.poly.ped"))


files2 <- matrix(paste0(rep(files, each = 2), c(".ped", ".map")),
                 ncol = 2, byrow = TRUE)

writeLines(paste(files2[, 1], files2[, 2])[-1],
           "../hapmap3_r2_b36_fwd.qc.poly/list_files.txt")



library(bigsnpr)

plink <- "../plink_linux_x86_64/plink"
bedfileQC <- snp_plinkQC(
  plink.path = plink,
  prefix.in = "../hapmap3_r2_b36_fwd.qc.poly/hapmap3",
  mind = 0.2,
  maf = 0.02,
  hwe = 1e-10,
  autosome.only = TRUE
)

bedfileQC2 <- snp_plinkIBDQC(
  plink.path = plink,
  bedfile.in = bedfileQC,
  pi.hat = 0.2,
  ncores = 3
)

map <- data.table::fread("../hapmap3_r2_b36_fwd.qc.poly/hapmap3_QC.bim",
                         data.table = FALSE)

celiac <- snp_attach("backingfiles/celiac300.rds")
mean(celiac$map$marker.ID %in% map$V2) # 78%
