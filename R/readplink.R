################################################################################

#' Read PLINK files into a "bigSNP"
#'
#' Function to read bed/bim/fam files into a [bigSNP][bigSNP-class].
#'
#' For more information on these formats, please visit
#' \href{http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#bed}{PLINK webpage}.
#' For other formats, please use PLINK to convert them in bedfiles,
#' which require minimal space to store and are faster to read.
#'
#' @param bedfile Path to file with extension ".bed" to read.
#' You need the corresponding ".bim" and ".fam" in the same directory.
#' @param backingfile The root name for the backing file(s) for the cache of
#' the [bigSNP][bigSNP-class] object.
#' @param backingpath The path to the directory containing the backing files.
#' Default is `"backingfiles"`.
#'
#' @return The path to one of the backing files which are created
#' by this function: \code{<backingfile>}.bk, \code{<backingfile>}.desc and
#' \code{<backingfile>}.rds in directory \code{<backingpath>}.\cr
#' __You shouldn't read from PLINK files more than once.__ Instead, use
#' [snp_attach] to load the "bigSNP" object in any R session from backing files.
#'
#' @note
#' The implementation was originally inspired from the code of
#' \href{https://github.com/andrewparkermorgan/argyle}{package argyle}.
#' I did many optimizations. Moreover, online reading
#' into a \code{big.matrix} makes it memory-efficient.
#'
#' @example examples/example.readplink.R
#' @export
snp_readBed <- function(bedfile,
                        block.size = 1000,
                        backingfile,
                        backingpath = "backingfiles") {
  backingpath <- path.expand(backingpath)
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
  names(fam) <- NAMES.FAM
  bim <- data.table::fread(bimfile, data.table = FALSE)
  names(bim) <- NAMES.MAP

  # prepare big.matrix
  n <- nrow(fam)
  m <- nrow(bim)
  bigGeno <- big.matrix(n, m, type = "char",
                        backingfile = paste0(backingfile, ".bk"),
                        backingpath = backingpath,
                        descriptorfile = paste0(backingfile, ".desc"))


  ## open bed file and check its magic number
  bed <- file(bedfile, open = "rb")
  magic <- readBin(bed, "raw", 3)
  if (!all(magic == c("6c", "1b", "01")))
    stop("Wrong magic number for bed file; should be -- 0x6c 0x1b 0x01 --.")

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

  rootPath <- file.path(backingpath, backingfile)
  saveRDS(snp_list, paste0(rootPath, ".rds"))

  paste0(rootPath, ".bk")
}

################################################################################

#' Attach a "bigSNP" from backing files
#'
#' Load a [bigSNP][bigSNP-class] from backing files into R.
#'
#' @param backingfile The path of one of the three (".bk", ".desc" or ".rds")
#' backing files for the cache of the "bigSNP" object.
#' @param readonly Is the \code{big.matrix} read only? Default is \code{TRUE}.
#'
#' @return The "bigSNP" object.
#' @example examples/example.readplink.R
#'
#' @export
snp_attach <- function(backingfile, readonly = TRUE) {
  root <- tools::file_path_sans_ext(backingfile)
  snp.list <- readRDS(paste0(root, ".rds"))

  snp.list$genotypes <- attach.big.matrix(paste0(root, ".desc"),
                                          readonly = readonly)

  snp.list$backingfile <- basename(root)
  snp.list$backingpath <- normalizePath(dirname(root))

  snp.list
}

################################################################################

#' Write PLINK files from a "bigSNP"
#'
#' Function to write bed/bim/fam files from a [bigSNP][bigSNP-class].
#'
#' @inheritParams bigsnpr-package
#' @param bedfile Path to file with extension ".bed" to create.
#'
#' @return The input `bedfile` path.
#'
#' @example
#' @export
snp_writeBed <- function(x, bedfile) {
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

#' Attach an "bigSNP" for examples and tests
#'
#' @inheritParams snp_readBed
#'
#' @return The example "bigSNP".
#'
#' @export
snp_attachExample <- function(backingfile = "test_doc",
                              backingpath = "backingfiles") {

  PATH <- file.path(backingpath, paste0(backingfile, ".bk"))
  if (!file.exists(PATH)) {
    bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
    PATH <- snp_readBed(bedfile, backingfile = backingfile,
                        backingpath = backingpath)
  }

  snp_attach(PATH)
}

################################################################################
