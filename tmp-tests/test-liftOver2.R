snp_modifyBuild <- function(info_snp, liftOver, from = "hg19", to = "hg18") {

  if (!all(c("chr", "pos") %in% names(info_snp)))
    stop("Please use proper names for variables in 'info_snp'. ",
         "Expected 'chr' and 'pos'.")

  # Need BED UCSC file for liftOver
  BED <- tempfile(fileext = ".BED")
  info_BED <- with(info_snp, data.frame(
    paste0("chr", chr), pos0 = pos - 1L, pos, id = rows_along(info_snp)))
  bigreadr::fwrite2(info_BED, BED, col.names = FALSE, sep = " ")

  # Make sure liftOver is executable
  Sys.chmod(liftOver, mode = (file.info(liftOver)$mode | "111"))

  # Need chain file
  url <- paste0("ftp://hgdownload.cse.ucsc.edu/goldenPath/", from, "/liftOver/",
                from, "To", tools::toTitleCase(to), ".over.chain.gz")
  chain <- tempfile(fileext = ".over.chain.gz")
  download.file(url, destfile = chain)

  # Run liftOver (usage: liftOver oldFile map.chain newFile unMapped)
  lifted <- tempfile(fileext = ".BED")
  unmaped <- tempfile(fileext = ".txt")
  system(paste(liftOver, BED, chain, lifted, unmaped))

  # readLines(lifted, n = 5)
  new_pos <- bigreadr::fread2(lifted)

  # readLines(unmaped, n = 6)
  bad <- grep("^#", readLines(unmaped), value = TRUE, invert = TRUE)
  message(length(bad), " variants have not been mapped.")

  new_info_snp <- info_snp[new_pos$V4, ]
  new_info_snp$pos <- new_pos$V3
  new_info_snp
}

bedfile.ref <- download_1000G(".")
bed.ref <- bed(bedfile.ref)
ref.map <- setNames(bed.ref$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
ref.map2 <- snp_modifyBuild(ref.map, "./liftOver")

bed.new <- bed("../../POPRES_data/POPRES_allchr.bed")
new.map <- setNames(bed.new$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
bigsnpr::snp_match(cbind(ref.map2, beta = 1), new.map, join_by_pos = FALSE)
nrow(.Last.value)  # 0
bigsnpr::snp_match(cbind(ref.map2, beta = 1), new.map)
nrow(.Last.value)  # 303,122

new.map2 <- snp_modifyBuild(new.map, "./liftOver", from = "hg18", to = "hg19")
bigsnpr::snp_match(cbind(ref.map, beta = 1), new.map2)
nrow(.Last.value)  # 304,953
