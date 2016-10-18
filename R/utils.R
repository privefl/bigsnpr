################################################################################

checkExists <- function(backingfile, backingpath) {
  # check directory and future file
  if (!file.exists(backingpath))
    stop(sprintf("Directory \"%s\" doesn't exist", backingpath))
  if (file.exists(file.path(backingpath, paste0(backingfile, ".bk"))))
    stop(sprintf("File \"%s.bk\" already exists in directory \"%s\"",
                 backingfile, backingpath))
}

################################################################################

check_x <- function(x, check.y = FALSE) {
  if (class(x) != "bigSNP") stop("x must be a bigSNP")

  X <- x$genotypes

  if (check.y) {
    y <- x$fam$pheno
    bigstatsr::check_X(X, y, y.type = "class")
  } else {
    bigstatsr::check_X(X)
  }
}

################################################################################

printf <- function(...) cat(sprintf(...))

################################################################################

CutBySize <- function(m, block.size) {
  nb <- ceiling(m / block.size)
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
