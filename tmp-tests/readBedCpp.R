Rcpp::sourceCpp('tmp-tests/fstream.cpp')
source('R/utils.R')
require(bigsnpr)

bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")

BedToBig2 <- function(bedfile,
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
  bigGeno <- bigmemory::big.matrix(n, m, type = "char",
                                   backingfile = paste0(backingfile, ".bk"),
                                   backingpath = backingpath,
                                   descriptorfile = paste0(backingfile, ".desc"))

  ## block size in bytes: (number of individuals)/4, to nearest byte
  bsz <- ceiling(n/4)

  ## open bed file and check its magic number
  bed <- file(bedfile, open = "rb")
  magic <- readBin(bed, "raw", 3)
  if (!all(magic == c("6c", "1b", "01"))) {
    stop("Wrong magic number for bed file; should be -- 0x6c 0x1b 0x01 --.")
  }
  close(bed)

  # read genotype data
  readbina(bedfile, bigGeno@address, bsz)

  snp_list <- list(genotypes = bigGeno, fam = fam, map = bim,
                   backingfile = backingfile,
                   backingpath = backingpath)
  class(snp_list) <- "bigSNP"

  saveRDS(snp_list, file.path(backingpath, paste0(backingfile, ".rds")))

  return(snp_list)
}

test <- BedToBig2(bedfile, "test3")
true <- AttachBigSNP("test_doc")

all.equal(test$genotypes[,], true$genotypes[,])


bedfile2 <- "../thesis-celiac/Dubois2010_data/FinnuncorrNLITUK3hap550.bed"
#"../stage-timc/Dubois2010_data/FinnuncorrNLITUK3hap550.bed"

print(system.time(
  test2 <- BedToBig2(bedfile2, backingfile = "test2")
))
# Windows: 10 min
# Linux: 30 sec

# print(system.time(
#   test3 <- BedToBig(bedfile2, block.size = 1e3, backingfile = "test7")
# ))
