wd <- setwd("~/Bureau/bigsnpr/tmp-data")

liftOverPlink <- "./liftOverPlink.py"
download.file("https://raw.githubusercontent.com/sritchie73/liftOverPlink/master/liftOverPlink.py",
              destfile = liftOverPlink)
Sys.chmod(liftOverPlink, mode = (file.info(liftOverPlink)$mode | "111"))
system(liftOverPlink)

from <- 19; to <- 18
chain <- paste0("./hg", from, "ToHg", to, ".over.chain.gz")
download.file(paste0("ftp://hgdownload.cse.ucsc.edu/goldenPath/hg", from,
                     "/liftOver/", basename(chain)),
              destfile = chain)

liftOver <- "./liftOver"
if (!file.exists(liftOver)) {
  download.file("http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver",
                destfile = liftOver)
  Sys.chmod(liftOver, mode = (file.info(liftOver)$mode | "111"))
}
system(liftOver)


bedfile.ref <- download_1000G(".")
bed.ref <- bed(bedfile.ref)
ref.map <- setNames(bed.ref$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
bigreadr::fwrite2(bed.ref$map[1:4], col.names = FALSE, sep = "\t",
                  file = (map <- tempfile(fileext = ".map")))

lifted <- "lifted"
system(glue::glue(
  "python {liftOverPlink} --map {map} --out {lifted} --chain {chain} --bin {liftOver}"))

system("python --version") # Python 2.7.5



library(dplyr)

# ref.map with positions in other build
ref.map2 <- bigreadr::fread2(paste0(lifted, ".map")) %>%
  mutate(V1 = as.integer(V1)) %>%
  left_join(ref.map, ., by = c(chr = "V1", rsid = "V2")) %>%
  select(chr, rsid, pos = V4, a0, a1)
str(ref.map2)

# POPRES with rsIDs
new.map2 <- bigreadr::fread2("../../POPRES_data/POPRES_Snps_QC2.txt", header = FALSE,
                 select = 1:2, col.names = c("id", "rsid")) %>%
  left_join(bigreadr::fread2("../../POPRES_data/POPRES_allchr.bim"), ., by = c(V2 = "id")) %>%
  select(chr = V1, rsid, pos = V4, a0 = V5, a1 = V6)
str(new.map2)

bigsnpr::snp_match(cbind(ref.map2, beta = 1), new.map2, join_by_pos = FALSE)

setwd(wd)
