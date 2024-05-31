################################################################################

# Same as in bigstatsr
CutBySize <- bigstatsr:::CutBySize
seq2      <- bigstatsr:::seq2
as_vec    <- bigstatsr:::as_vec

################################################################################

# global variable definitions due to non standard evaluations
utils::globalVariables(c("ic", "f", "lp", "LD.wiki34", "OS", "arch", "avx2",
                         "chromosome", "position", "allele1", "allele2",
                         "marker.ID", "rsid", "chr", "pos", "p_init",
                         "df_af", ".", "FREQ", "ID", "INFO", "ind_rm",
                         "file_start_position", "h2", "p", "sparse",
                         "arm", "amd", "id_current"))

################################################################################

# functions for encoding/decoding bed files
getCode <- function(NA.VAL = 3L) {
  geno.raw <- !as.logical(rawToBits(as.raw(0:255)))
  s <- c(TRUE, FALSE)
  geno1 <- geno.raw[s]
  geno2 <- geno.raw[!s]
  geno <- geno1 + geno2
  geno[!geno1 & geno2] <- NA.VAL
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

write.table2 <- function(x, file, ...) {
  data.table::fwrite(x, file, sep = "\t", quote = FALSE, scipen = 50,
                     row.names = FALSE, col.names = FALSE, ...)
}

################################################################################

getNewFile <- function(x, type) {

  root <- sub_bk(x$genotypes$backingfile)
  EXTS <- c("bk", "rds")

  number <- 1
  repeat {
    files <- sprintf("%s_%s%d.%s", root, type, number, EXTS)
    if (all(!file.exists(files))) break
    number <- number + 1
  }

  sprintf("%s_%s%d", root, type, number)
}

################################################################################
