################################################################################

# functions for encoding/decoding bed files
getCode <- function(NA.VAL = 3L) {
  all.raws <- as.raw(0:255)
  geno.raw <- as.logical(rawToBits(all.raws))
  s <- c(TRUE, FALSE)
  geno1 <- geno.raw[s]
  geno2 <- geno.raw[!s]
  geno <- geno1 + geno2
  geno[geno1 & !geno2] <- NA.VAL
  dim(geno) <- c(4, 256)
  storage.mode(geno) <- "raw"
  geno
}
# t(mapply(rep, times = 4^(3:0), each = 4^(0:3),
#          MoreArgs = list(x = as.raw(c(0, 3, 1, 2)))))

getInverseCode <- function() {
  geno <- getCode()
  storage.mode(geno) <- "integer"
  r <- raw(256)
  dim(r) <- rep(4, 4)
  for (i in 1:256) {
    ind <- geno[, i] + 1
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

printf <- function(...) cat(sprintf(...))

message2 <- function(...) message(sprintf(...))

stop2 <- function(...) stop(sprintf(...), call. = FALSE)

################################################################################

write.table2 <- function(x, file) {
  utils::write.table(x = x, file = file, quote = FALSE, sep = "\t",
                     row.names = FALSE, col.names = FALSE)
}

################################################################################

dir.create2 <- function(dir.path) {
  if (!dir.exists(dir.path)) {
    if (dir.create(dir.path)) {
      message2("Creating directory \"%s\" which didn't exist", dir.path)
    } else {
      stop2("Problem creating directory \"%s\". Recursive path?", dir.path)
    }
  }
}

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
