source('R/utils.R')

BedToBig4 <- function(bedfile,
                      block.size,
                      backingfile,
                      backingpath = "backingfiles") {
  checkExists(backingfile, backingpath)

  # check extension of file
  ext <- tools::file_ext(bedfile)
  if (ext != "bed") {
    stop(sprintf("Extension .%s unsupported, requires .bed instead", ext))
  } else {
    bimfile <- sub("\\.bed$", ".bim", bedfile)
    famfile <- sub("\\.bed$", ".fam", bedfile)
  }

  # check if all three files exist
  if (!file.exists(bedfile)) {
    stop(sprintf("File \"%s\" doesn't exist", bedfile))
  } else if (!file.exists(bimfile)) {
    stop(sprintf("File \"%s\" doesn't exist", bimfile))
  } else if (!file.exists(famfile)) {
    stop(sprintf("File \"%s\" doesn't exist", famfile))
  }

  # read map and family files
  fam <- data.table::fread(famfile, data.table = FALSE)
  names(fam) <- c("family.ID", "sample.ID", "paternal.ID",
                  "maternal.ID", "sex", "affection")
  bim <- data.table::fread(bimfile, data.table = FALSE)
  names(bim) <- c("chromosome", "marker.ID", "genetic.dist",
                  "physical.pos", "allele1", "allele2")

  # prepare big.matrix
  n <- nrow(fam)
  m <- nrow(bim)
  bsz <- ceiling(n/4)
  bigGeno <- bigmemory::big.matrix(bsz, m, type = "raw",
                                   backingfile = paste0(backingfile, ".bk"),
                                   backingpath = backingpath,
                                   descriptorfile = paste0(backingfile, ".desc"))


  ## open bed file and check its magic number
  bed <- file(bedfile, open = "rb")
  magic <- readBin(bed, "raw", 3)
  if (!all(magic == c("6c", "1b", "01"))) {
    stop("Wrong magic number for bed file; should be -- 0x6c 0x1b 0x01 --.")
  }

  # now actually read genotypes block by block
  intervals <- CutBySize(m, block.size)
  nb.blocks <- nrow(intervals)
  if (intr <- interactive())
    pb <- utils::txtProgressBar(min = 0, max = nb.blocks, style = 3)

  colOffset <- 0L
  for (k in 1:nb.blocks) {
    size <- intervals[k, "size"]

    bigGeno[, 1:size + colOffset] <-
      readBin(bed, what = "raw", n = bsz * size)

    colOffset <- colOffset + size
    if (intr) utils::setTxtProgressBar(pb, k)
  }
  close(bed)
  if (intr) close(pb)


  snp_list <- list(genotypes = bigGeno, fam = fam, map = bim,
                   backingfile = backingfile,
                   backingpath = backingpath)
  class(snp_list) <- "bigSNP"

  saveRDS(snp_list, file.path(backingpath, paste0(backingfile, ".rds")))

  snp_list
}

bedfile2 <- "../stage-timc/Dubois2010_data/FinnuncorrNLITUK3hap550.bed"

print(system.time(
  test <- BedToBig4(bedfile2, block.size = 500, backingfile = "test6")
))
# 74 sec just to read binary to logical
# 89 sec write raw data
