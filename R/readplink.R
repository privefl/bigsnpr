################################################################################

#'@title Read PLINK files into a "bigSNP".
#'@description Functions to read ped/map or bed/bim/fam files
#'into a \code{bigSNP}.\cr
#'For more information on these formats, please visit
#'\href{http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped}{PLINK webpage}.
#'Prefer using bedfiles than pedfiles because
#' they require minimal space to store and are faster to read.
#'@param bedfile Path to file with extension .bed. You need the corresponding
#'.bim and .fam in the same directory.
#'@param pedfile Path to file with extension .ped. You need the corresponding
#'.map in the same directory.
#'@param block.size Maximum number of loci read at once (for all individuals).
#'@param backingfile The root name for the backing file(s) for the cache of
#'the resulting object.
#'@param backingpath The path to the directory containing the file backing cache.
#'Default is "backingfiles". It needs to exist.
#'@param readonly Is the \code{big.matrix} read only ? Default is \code{TRUE}.
#'@return A \code{bigSNP}.\cr
#'Reading PLINK files creates
#'\code{backingfile}.bk, \code{backingfile}.desc and \code{backingfile}.rds
#'in directory \code{backingpath}.\cr
#'You shouldn't read from PLINK files more than once.
#'Instead, use \code{AttachBigSNP}
#'to load this object in another session from backing files.
#'@note
#'The implementation is inspired from the code of
#'\href{https://github.com/andrewparkermorgan/argyle}{package argyle},
#'with some optimizations.\cr
#'Especially, online reading into a \code{big.matrix} makes it memory-efficient.
#'@examples \dontrun{
#'
#'bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
#'
#'if (!dir.exists("backingfiles"))
#'  dir.create("backingfiles")
#'if (file.exists("backingfiles/test_doc"))
#'  file.remove("backingfiles/test_doc")
#'if (file.exists("backingfiles/test_doc.desc"))
#'  file.remove("backingfiles/test_doc.desc")
#'if (file.exists("backingfiles/test_doc.rds"))
#'  file.remove("backingfiles/test_doc.rds")
#'
#'# Reading the bedfile and storing the data in directory "backingfiles"
#'test <- BedToBig(bedfile, 50, "test_doc")
#'
#'# Removing the R object
#'rm(test)
#'
#'# Loading it from backing files
#'test <- AttachBigSNP("test_doc")
#'
#'str(test)
#'print(dim(test$genotypes))
#'print(test$genotypes[1:8, 1:8])
#'}
#'@seealso \code{\link{bigSNP}} \code{\link{dir.create}}
#'@name readplink
NULL


################################################################################

#' @rdname readplink
#' @export
BedToBig <- function(bedfile,
                     block.size,
                     backingfile,
                     backingpath = "backingfiles") {
  if (!file.exists(backingpath))
    stop(sprintf("Directory \"%s\" doesn't exist", backingpath))
  if ( file.exists(file.path(backingpath, backingfile)))
    stop(sprintf("File \"%s\" already exists in directory \"%s\"",
                 backingfile, backingpath))

  ListToInd <- function(list, colOffset) {
    cbind(row = unlist(list),
          col = rep(seq_along(list), sapply(list, length)) + colOffset)
  }

  # check extension of file
  ext <- tools::file_ext(bedfile)
  if (ext != "bed") {
    stop(sprintf("Extension .%s unsupported, requires .bed instead", ext))
  } else {
    bimfile <- sub(".bed$", ".bim", bedfile)
    famfile <- sub(".bed$", ".fam", bedfile)
  }

  # check if all three files exist
  if (!file.exists(bedfile)) {
    stop(sprintf("File \"%s\" doesn't exist\n", bedfile))
  } else if (!file.exists(bimfile)) {
    stop(sprintf("File \"%s\" doesn't exist\n", bimfile))
  } else if (!file.exists(famfile)) {
    stop(sprintf("File \"%s\" doesn't exist\n", famfile))
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

  ## now actually read genotypes block by block
  intervals <- CutBySize(m, block.size)
  nb.blocks <- nrow(intervals)
  if (intr <- interactive())
    pb <- utils::txtProgressBar(min = 0, max = nb.blocks, style = 3)
  s1 <- seq(1, 2*n, 2)
  s2 <- s1 + 1

  colOffset <- 0
  for (k in 1:nb.blocks) {
    if (intr) utils::setTxtProgressBar(pb, k - 1)
    list.ind.na <- list()
    size <- intervals[k, "size"]
    geno.mat <- matrix(0, n, size)
    for (i in 1:size) {
      geno.raw <- as.logical(rawToBits(readBin(bed, "raw", bsz)))
      geno1 <- geno.raw[s1]
      geno2 <- geno.raw[s2]
      ## express genotypes as minor allele dosage (0,1,2)
      geno.mat[, i] <- geno1 + geno2
      ## recall that 0/1 is het, but 1/0 is missing
      list.ind.na[[i]] <- which(geno1 & !geno2)
    }
    rawToBigPart(geno.mat, bigGeno@address, colOffset)
    ind.na <- ListToInd(list.ind.na, colOffset)
    if (nrow(ind.na) > 0) bigGeno[ind.na] <- NA
    colOffset <- colOffset + size
  }
  close(bed)
  if (intr) {
    utils::setTxtProgressBar(pb, nb.blocks)
    close(pb)
  }

  snp_list <- list(genotypes = bigGeno, fam = fam, map = bim,
                   backingfile = backingfile,
                   backingpath = backingpath)
  class(snp_list) <- "bigSNP"

  saveRDS(snp_list, file.path(backingpath, paste0(backingfile, ".rds")))

  return(snp_list)
}

################################################################################

#' @rdname readplink
#' @export
PedToBig <- function(pedfile,
                     block.size,
                     backingfile,
                     backingpath = "backingfiles") {
  if (!file.exists(backingpath))
    stop(sprintf("Directory \"%s\" doesn't exist", backingpath))
  if ( file.exists(file.path(backingpath, backingfile)))
    stop(sprintf("File \"%s\" already exists in directory \"%s\"",
                 backingfile, backingpath))

  dna.letters <- c("A", "C", "T", "G")

  CompareToRef <- function(x, ref, s) {
    comp1 <- ifelse(x[s]     %in% dna.letters, (x[s] != ref), NA)
    comp2 <- ifelse(x[s + 1] %in% dna.letters, (x[s + 1] != ref), NA)

    return(comp1 + comp2)
  }

  diffRef <- function(letter, tmp, s) {
    rowSums(sapply(tmp, CompareToRef, ref = letter, s = s), na.rm = T)
  }

  # check extension of file
  ext <- tools::file_ext(pedfile)
  if (ext != "ped") {
    stop(sprintf("Extension .%s unsupported, requires .ped instead", ext))
  } else {
    mapfile <- sub(".ped$", ".map", pedfile)
  }

  # check if both files exist
  if (!file.exists(pedfile)) {
    stop(sprintf("File \"%s\" doesn't exist\n", pedfile))
  } else if (!file.exists(mapfile)) {
    stop(sprintf("File \"%s\" doesn't exist\n", mapfile))
  }

  # get the number of SNPs
  ped <- file(pedfile, open = "r")
  tmp <- readLines(ped, n = 1)
  close(ped)
  tmp <- strsplit(tmp, " ", fixed = T)
  m.all <- length(tmp[[1]])
  s <- seq(7L, m.all, 2L)
  m <- length(s)
  map <- data.table::fread(mapfile, data.table = FALSE)
  if (m != nrow(map)) {
    stop(sprintf("%d markers were read from the .map file and
                 %d from the .ped file", nrow(map), m))
  }
  #data.table::setnames(map, 1:4, c("chromosome", "marker.ID",
  #                                 "genetic.dist", "physical.pos"))
  names(map) <- c("chromosome", "marker.ID", "genetic.dist",
                  "physical.pos")

  # first read to get useful infos
  printf("\nBegin first read to get useful info\n")
  counts <- matrix(0L, m, 4)
  n <- 0L
  ped <- file(pedfile, open = "r")
  while ((n.part <- length(tmp <- readLines(ped, n = block.size))) > 0) {
    tmp <- strsplit(tmp, " ", fixed = T)
    counts <- counts + sapply(dna.letters, diffRef, tmp = tmp, s = s)
    n <- n + n.part
    printf("%d lines have been read so far\n", n)
  }
  close(ped)
  ind.rm <- which(counts == 0, arr.ind = T)[, 1]
  if ((len <- length(ind.rm)) > 0) {
    printf("\nRemoving %d monoallelic markers\n", len)
    map <- map[-ind.rm]
    alleles <- apply(counts[-ind.rm, ], 1, order)
  } else {
    alleles <- apply(counts, 1, order)
  }
  ref <- dna.letters[apply(counts, 1, which.min)]
  rm(counts)
  map$allele1 <- dna.letters[alleles[1, ]]
  map$allele2 <- dna.letters[alleles[2, ]]
  rm(alleles)

  # second read to fill the genotypic matrix
  intervals <- CutBySize(n, block.size)
  nb.blocks <- nrow(intervals)
  printf("\nSecond and last read to fill the genotypic matrix\n")
  if (intr <- interactive())
    pb <- utils::txtProgressBar(min = 0, max = nb.blocks + 1, style = 3)
  bigGeno <- bigmemory::big.matrix(n, m - len, type = "char",
                                   backingfile = paste0(backingfile, ".bk"),
                                   backingpath = backingpath,
                                   descriptorfile = paste0(backingfile, ".desc"))
  fam <- list()
  opt.save <- options(bigmemory.typecast.warning = FALSE)
  ped <- file(pedfile, open = "r")
  for (k in 1:nb.blocks) {
    if (intr) utils::setTxtProgressBar(pb, k - 1)
    tmp <- readLines(ped, n = intervals[k, "size"])
    tmp <- strsplit(tmp, " ", fixed = T)
    fam[[k]] <-  sapply(tmp, head, n = 6L)
    if (len > 0) {
      tmp <- t(sapply(tmp, CompareToRef, ref = ref, s = s)[-ind.rm, ])
    } else {
      tmp <- t(sapply(tmp, CompareToRef, ref = ref, s = s))
    }
    bigGeno[seq2(intervals[k, ]), ] <- tmp
  }
  close(ped)
  options(opt.save)
  if (intr) utils::setTxtProgressBar(pb, nb.blocks)

  # convert the fam list in a data.frame
  tmpfile.name <- tempfile()
  cat("", file = tmpfile.name)
  lapply(fam, function(l) apply(l, 2, function(text)
    cat(text, file = tmpfile.name, fill = TRUE, append = TRUE)))
  fam <- data.table::fread(tmpfile.name, data.table = FALSE)
  names(fam) <- c("family.ID", "sample.ID", "paternal.ID",
                  "maternal.ID", "sex", "affection")
  #file.remove(file = tmpfile.name)

  if (intr) {
    utils::setTxtProgressBar(pb, nb.blocks + 1)
    close(pb)
  }

  snp_list <- list(genotypes = bigGeno, fam = fam, map = map,
                   backingfile = backingfile,
                   backingpath = backingpath)
  class(snp_list) <- "bigSNP"

  saveRDS(snp_list, file.path(backingpath, paste0(backingfile, ".rds")))

  return(snp_list)
}

################################################################################

#' @rdname readplink
#' @export
AttachBigSNP <- function(backingfile,
                         backingpath = "backingfiles",
                         readonly = TRUE) {
  snp.list <- readRDS(file.path(backingpath, paste0(backingfile, ".rds")))

  snp.list$genotypes <-
    bigmemory::attach.big.matrix(file.path(backingpath,
                                           paste0(backingfile, ".desc")),
                                 readonly = readonly)

  snp.list$backingfile <- backingfile
  snp.list$backingpath <- backingpath

  return(snp.list)
}

################################################################################
