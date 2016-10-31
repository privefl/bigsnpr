#Rcpp::sourceCpp('src/readplink.cpp')
source('R/utils.R')

require(bigmemory)

BedToBig5 <- function(bedfile,
                      block.size,
                      backingfile,
                      backingpath = "backingfiles", ncores) {
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
  X.desc <- describe(bigGeno)

  ## open bed file and check its magic number
  bed <- file(bedfile, open = "rb")
  magic <- readBin(bed, "raw", 3)
  if (!all(magic == c("6c", "1b", "01"))) {
    stop("Wrong magic number for bed file; should be -- 0x6c 0x1b 0x01 --.")
  }
  close(bed)

  # match each possible code
  getCode <- function() {
    all.raws <- as.raw(0:255)
    geno.raw <- as.logical(rawToBits(all.raws))
    s <- c(TRUE, FALSE)
    geno1 <- geno.raw[s]
    geno2 <- geno.raw[!s]
    geno <- geno1 + geno2
    geno[geno1 & !geno2] <- NA_CHAR
    dim(geno) <- c(4, 256)
    geno
  }
  geno <- getCode()

  ## block size in bytes: (number of individuals)/4, to nearest byte
  bsz <- ceiling(n/4)
  intervals <- CutBySize(m, block.size)

  part <- function(i) {
    X <- attach.big.matrix(X.desc, backingpath = backingpath)

    # reading
    bed <- file(bedfile, open = "rb")
    magic <- readBin(bed, "raw", 3)

    # now actually read genotypes block by block
    nb.blocks <- nrow(intervals)
    colOffset <- 0
    for (k in 1:nb.blocks) {
      size <- intervals[k, "size"]
      # everyone read
      tmp <- readBin(bed, "raw", bsz * size)

      # determine which will be writing this
      which.core <- (k %% ncores) + 1
      if (which.core == i) {
        rawToBigPart(X@address,
                     source = tmp,
                     tab = geno,
                     size = size, colOffset = colOffset,
                     n = n, bsz = bsz)
      }
      colOffset <- colOffset + size
    }

    close(bed)
    0
  }

  obj <- foreach::foreach(i = 1:ncores, .packages = "bigsnpr")
  expr_fun <- function(i) part(i)
  foreach2(obj, expr_fun, ncores)



  snp_list <- list(genotypes = bigGeno, fam = fam, map = bim,
                   backingfile = backingfile,
                   backingpath = backingpath)
  class(snp_list) <- "bigSNP"

  saveRDS(snp_list, file.path(backingpath, paste0(backingfile, ".rds")))

  snp_list
}

bedfile2 <- "../stage-timc/Dubois2010_data/FinnuncorrNLITUK3hap550.bed"

print(system.time(
  test <- BedToBig5(bedfile2, block.size = 500, backingfile = "test2",
                    ncores = 4)
))
# 74 sec just to read binary to logical
# 308 with 4 cores...
