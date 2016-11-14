################################################################################

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
  if (!dir.exists(backingpath))
    stop(sprintf("Directory \"%s\" doesn't exist", backingpath))
  if (file.exists(file.path(backingpath, paste0(backingfile, ".bk"))))
    stop(sprintf("File \"%s.bk\" already exists in directory \"%s\"",
                 backingfile, backingpath))
}

################################################################################

check_x <- function(x, check.y = FALSE) {
  if (class(x) != "bigSNP") stop("x must be a bigSNP.")

  X <- x$genotypes
  if (class(X) != "big.matrix") stop("X must be a big.matrix.")

  if (check.y) {
    y <- x$fam$pheno
    if (is.null(y))
      stop("Please use _GetPhenos_ before using this function.")
    if (!isTRUE(all.equal(sort(unique(y)), c(-1, 1))))
      stop("y should be a vector of 1 (cases) and -1 (controls).")
  }
}

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

foreach2 <- function(obj, expr_fun, ncores, outfile = NULL) {
  if (is.seq <- (ncores == 1)) {
    foreach::registerDoSEQ()
  } else {
    if (is.null(outfile)) {
      cl <- parallel::makeCluster(ncores)
    } else {
      cl <- parallel::makeCluster(ncores, outfile = outfile)
    }
    doParallel::registerDoParallel(cl)
  }
  res <- eval(parse(text = sprintf("foreach::`%%dopar%%`(obj, expr_fun(%s))",
                                   obj$argnames)))
  if (!is.seq) parallel::stopCluster(cl)

  return(res)
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
