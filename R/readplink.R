################################################################################

#' @title Read PLINK files into a "bigSNP".
#' @description Functions to read bed/bim/fam files
#' into a `bigSNP`.\cr
#' For more information on these formats, please visit
#' \href{http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#bed}{PLINK webpage}.
#' For other formats, please use PLINK to convert them in bedfiles,
#' which require minimal space to store and are faster to read.
#' @param bedfile Path to file with extension .bed. You need the corresponding
#' .bim and .fam in the same directory.
#' @param backingfile The root name for the backing file(s) for the cache of
#' the resulting object.
#' @param backingpath The path to the directory containing the file backing cache.
#' Default is "backingfiles". It needs to exist (use [dir.create]).
#' @param readonly Is the \code{big.matrix} read only? Default is \code{TRUE}.
#' @return A \code{bigSNP}.\cr
#' Reading PLINK files creates
#' \code{backingfile}.bk, \code{backingfile}.desc and \code{backingfile}.rds
#' in directory \code{backingpath}.\cr
#' You shouldn't read from PLINK files more than once.
#' Instead, use \code{AttachBigSNP}
#' to load this object in another session from backing files.
#' @note
#' The implementation was originally inspired from the code of
#' \href{https://github.com/andrewparkermorgan/argyle}{package argyle}.
#' I did many optimizations. Especially, online reading
#' into a \code{big.matrix} makes it memory-efficient.
#' @example examples/example.readplink.R
#' @name readplink
NULL

################################################################################

#' @rdname readplink
#' @export
BedToBig <- function(bedfile,
                     block.size = 1000,
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


  ## open bed file and check its magic number
  bed <- file(bedfile, open = "rb")
  magic <- readBin(bed, "raw", 3)
  if (!all(magic == c("6c", "1b", "01"))) {
    stop("Wrong magic number for bed file; should be -- 0x6c 0x1b 0x01 --.")
  }

  # match each possible code
  geno <- getCode()

  ## block size in bytes: (number of individuals)/4, to nearest byte
  bsz <- ceiling(n/4)


  # now actually read genotypes block by block
  intervals <- CutBySize(m, block.size)
  nb.blocks <- nrow(intervals)

  colOffset <- 0
  for (k in 1:nb.blocks) {
    size <- intervals[k, "size"]
    rawToBigPart(bigGeno@address,
                 source = readBin(bed, "raw", bsz * size),
                 tab = geno,
                 size = size, colOffset = colOffset,
                 n = n, bsz = bsz)
    colOffset <- colOffset + size
  }
  close(bed)


  snp_list <- list(genotypes = bigGeno, fam = fam, map = bim,
                   backingfile = backingfile,
                   backingpath = normalizePath(backingpath))
  class(snp_list) <- "bigSNP"

  saveRDS(snp_list, file.path(backingpath, paste0(backingfile, ".rds")))

  snp_list
}

################################################################################

#' @rdname readplink
#' @param backingfile The path of one of the three (.bk, .desc or .rds)
#' backing files for the cache of the resulting object.
#' @export
AttachBigSNP <- function(backingfile, readonly = TRUE) {
  root <- tools::file_path_sans_ext(backingfile)
  snp.list <- readRDS(paste0(root, ".rds"))

  snp.list$genotypes <- attach.big.matrix(paste0(root, ".desc"),
                                          readonly = readonly)

  snp.list$backingfile <- basename(root)
  snp.list$backingpath <- normalizePath(dirname(root))

  return(snp.list)
}

################################################################################

#' Title
#'
#' @param x
#' @param bedfile
#'
#' @return
#' @export
#'
#' @examples
BigToBed <- function(x, bedfile) {
  check_x(x)

  # check extension of file
  ext <- tools::file_ext(bedfile)
  if (ext != "bed") {
    stop(sprintf("Extension .%s unsupported, requires .bed instead", ext))
  } else {
    bimfile <- sub("\\.bed$", ".bim", bedfile)
    famfile <- sub("\\.bed$", ".fam", bedfile)
  }

  if (file.exists(bedfile))
    stop(sprintf("File %s already exists", bedfile))

  # write map and family files
  write.table(x$fam[NAMES.FAM], file = famfile, quote = FALSE,
              sep = "\t", row.names = FALSE, col.names = FALSE)
  write.table(x$map[NAMES.MAP], file = bimfile, quote = FALSE,
              sep = "\t", row.names = FALSE, col.names = FALSE)

  ## write bed file
  X <- x$genotypes
  writebina(bedfile, X@address, getInverseCode())

  bedfile
}

################################################################################

#' @rdname readplink
#'
#' @export
snp_readExample <- function(backingfile = "test_doc",
                            backingpath = "backingfiles") {
  # Creating directory for backing files
  dir.create2(backingpath)

  path <- file.path(backingpath, paste0(backingfile, ".bk"))
  unlink2(path) # delete if exists

  # Reading the bedfile and storing the data
  bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
  test <- BedToBig(bedfile, backingfile = backingfile,
                   backingpath = backingpath)
}

################################################################################

#' @rdname readplink
#' @export
snp_readHapMap3 <- function(backingfile = "test_HapMap3",
                            backingpath = "backingfiles") {
  # Creating directory for backing files
  dir.create2(backingpath)

  path <- file.path(backingpath, paste0(backingfile, ".bk"))
  unlink2(path) # delete if exists

  # Reading the bedfile and storing the data
  bedfile <- system.file("extdata", "HapMap3.bed", package = "bigsnpr")
  test <- BedToBig(bedfile, backingfile = backingfile,
                   backingpath = backingpath)
}

################################################################################
