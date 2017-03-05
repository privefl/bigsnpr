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
                        backingpath = "backingfiles",
                        cpp = FALSE) {
  backingpath <- path.expand(backingpath)
  rootPath <- file.path(backingpath, backingfile)

  # check if file and path exist
  if (dir.exists(backingpath)) {
    if (file.exists(bkfile <- paste0(rootPath, ".bk"))) {
      message(sprintf("File \"%s.bk\" already exists in directory \"%s\"",
                      backingfile, backingpath))
      message("Aborting and returning its path..")
      return(bkfile)
    }
  } else {
    if (dir.create(backingpath))
      message(sprintf("Creating directory \"%s\" which didn't exist",
                      backingpath))
  }

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
  bigGeno <- big.matrix(n, m, type = "raw",
                        backingfile = paste0(backingfile, ".bk"),
                        backingpath = backingpath,
                        descriptorfile = paste0(backingfile, ".desc"))


  if (cpp) {
    readbina(bedfile, bigGeno@address, getCode())
  } else {
    ## open bed file and check its magic number
    bed <- file(bedfile, open = "rb")
    magic <- readBin(bed, "raw", 3)
    if (!all(magic == c("6c", "1b", "01")))
      stop("Wrong magic number for bed file; should be -- 0x6c 0x1b 0x01 --.")

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
                   tab = getCode(),
                   size = size, colOffset = colOffset,
                   n = n, bsz = bsz)
      colOffset <- colOffset + size
    }
    close(bed)
  }

  code <- rep(NA_real_, 256)
  code[1:3] <- c(0, 1, 2)
  bigGeno.code <- as.BM.code(bigGeno, code)

  rds <- paste0(rootPath, ".rds")

  snp_list <- structure(list(genotypes = describe(bigGeno.code),
                             fam = fam,
                             map = bim,
                             savedIn = paste0(rootPath, ".rds")),
                        class = "bigSNP")

  saveRDS(snp_list, rds)

  rds
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
#' @example examples/example-writeplink.R
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
  X <- attach.BM(x$genotypes)
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
snp_attachExtdata <- function(backingfile = "test_doc",
                              backingpath = "backingfiles") {

  PATH <- file.path(backingpath, paste0(backingfile, ".bk"))
  if (!file.exists(PATH)) {
    bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
    PATH <- snp_readBed(bedfile, backingfile = backingfile,
                        backingpath = backingpath)
  }

  readRDS(sub("\\.bk$", ".rds", PATH))
}

################################################################################
