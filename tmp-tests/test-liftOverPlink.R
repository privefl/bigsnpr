download.file("https://raw.githubusercontent.com/sritchie73/liftOverPlink/master/liftOverPlink.py",
              destfile = (liftOverPlink <- tempfile(fileext = ".py")))
Sys.chmod(liftOverPlink, mode = (file.info(liftOverPlink)$mode | "111"))
system(liftOverPlink)

download.file("ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
              destfile = (chain <- tempfile(fileext = ".over.chain.gz")))

download.file("http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver",
              destfile = (liftOver <- tempfile()))
Sys.chmod(liftOver, mode = (file.info(liftOver)$mode | "111"))
system(liftOver)

bim <- "../POPRES_data/POPRES_allchr.bim"
# download.file("https://raw.githubusercontent.com/gabraham/flashpca/master/HapMap3/data.bim",
#               destfile = (bim <- tempfile(fileext = ".bim")))
bigreadr::fwrite2(bigreadr::fread2(bim, select = 1:4), col.names = FALSE, sep = "\t",
                  file = (map <- tempfile(fileext = ".map")))

lifted <- tempfile()
system(glue::glue(
  "python {liftOverPlink} --map {map} --out {lifted} --chain {chain} --bin {liftOver}"))

system("python --version") # Python 2.7.5



bedfile.ref <- download_1000G("tmp-data")
bed.ref <- bed(bedfile.ref)
ref.map <- setNames(bed.ref$map[-3], c("chr", "rsid", "pos", "a1", "a0"))


map <- bigreadr::fread2(paste0(lifted, ".map"))
bim <- bigreadr::fread2(bim)
library(dplyr)
left_join(bim, mutate(map, V1 = as.integer(V1)), by = c("V1", "V2")) %>%
  select(chr = V1, pos = V4.y, a0 = V5, a1 = V6) %>%
  bigsnpr::snp_match(cbind(ref.map, beta = 1), .)
