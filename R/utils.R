################################################################################

# functions for encoding/decoding bed files
getCode <- function(NA_CHAR = -128L) {
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

getInverseCode <- function() {
  geno <- getCode(3) + 1
  r <- raw(256)
  dim(r) <- rep(4, 4)
  for (i in 1:256) {
    ind <- geno[, i]
    r[ind[1], ind[2], ind[3], ind[4]] <- as.raw(i - 1)
  }
  r
}

################################################################################

NAMES.MAP <- c("chromosome", "marker.ID", "genetic.dist",
               "physical.pos", "allele1", "allele2")

NAMES.FAM <- c("family.ID", "sample.ID", "paternal.ID",
               "maternal.ID", "sex", "affection")

################################################################################

checkExists <- function(backingfile, backingpath) {
  # check directory and future file
  if (dir.exists(backingpath)) {
    if (file.exists(file.path(backingpath, paste0(backingfile, ".bk"))))
      stop(sprintf("File \"%s.bk\" already exists in directory \"%s\"",
                   backingfile, backingpath))
  } else {
    if (dir.create(backingpath))
      message(sprintf("Creating directory \"%s\" which didn't exist",
                      backingpath))
  }
}

################################################################################

check_x <- function(x) {
  if (class(x) != "bigSNP") stop("x must be a bigSNP.")
  if (class(x$genotypes) != "big.matrix")
    stop("x$genotypes must be a big.matrix.")
}

################################################################################

# dir.create2 <- function(path) {
#   if (!dir.exists(path)) {
#     printf("Creating directory \"%s\"..\n", path)
#     dir.create(path)
#   }
# }
#
# unlink2 <- function(path) {
#   if (file.exists(path)) unlink(path)
# }

################################################################################

printf <- function(...) cat(sprintf(...))

################################################################################

CutBySize <- function(m, block.size, nb = ceiling(m / block.size)) {
  int <- m / nb

  upper <- round(1:nb * int)
  lower <- c(1, upper[-nb] + 1)
  size <- c(upper[1], diff(upper))

  cbind(lower, upper, size)
}

################################################################################

seq2 <- function(lims) {
  seq(lims[1], lims[2])
}

################################################################################

LimsChr <- function(infos) {
  map.rle <- rle(infos$map$chromosome)
  upper <- cumsum(map.rle$length)
  lower <- c(1, upper[-length(upper)] + 1)

  cbind(lower, upper, "chr" = map.rle$values)
}

################################################################################

checkFile <- function(x, type) {
  number <- 1
  while (file.exists(file.path(
    x$backingpath,
    paste0(newfile <- paste0(x$backingfile, "_", type, number), ".bk")))) {
    number <- number + 1
  }

  newfile
}

################################################################################

transform_levels <- function(y, new.levels = 0:1) {
  y2 <- factor(y, ordered = TRUE)
  lvl <- levels(y2)
  if (length(lvl) != 2)
    stop("You must have exactly two levels in y.")
  levels(y2) <- new.levels
  as.numeric(as.character(y2))
}

################################################################################
