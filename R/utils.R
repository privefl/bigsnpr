################################################################################
# Same as in bigstatsr

printf <- bigstatsr:::printf

message2 <- bigstatsr:::message2

stop2 <- bigstatsr:::stop2

transform_levels <- bigstatsr:::transform_levels

CutBySize <- bigstatsr:::CutBySize

seq2 <- bigstatsr:::seq2

################################################################################

# functions for encoding/decoding bed files
getCode <- function(NA.VAL = 3L) {
  geno.raw <- as.logical(rawToBits(RAWS))
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

write.table2 <- function(x, file) {
  utils::write.table(x = x, file = file,
                     quote = FALSE,
                     sep = "\t",
                     row.names = FALSE,
                     col.names = FALSE)
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

LimsChr <- function(infos) {
  map.rle <- rle(infos$map$chromosome)
  upper <- cumsum(map.rle$length)
  lower <- c(1, upper[-length(upper)] + 1)

  cbind(lower, upper, "chr" = map.rle$values)
}

################################################################################

getNewFiles <- function(rdsfile, type) {

  root <- sub("\\.rds$", "", rdsfile)
  EXTS <- c("bk", "desc", "rds")

  number <- 1
  repeat {
    files <- sprintf("%s_%s%d.%s", root, type, number, EXTS)
    if (all(!file.exists(files))) break
    number <- number + 1
  }

  as.list(structure(files, names = EXTS))
}

################################################################################
