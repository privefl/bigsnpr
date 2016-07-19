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

foreach2 <- function(obj, expr_fun, ncores) {
  if (is.seq <- (ncores == 1)) {
    foreach::registerDoSEQ()
  } else {
    cl <- parallel::makeCluster(ncores)
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
  while (file.exists(
    file.path(x$backingpath,
              paste0(newfile <- paste0(x$backingfile, "_", type, number),
                     ".desc")))) {
    number <- number + 1
  }

  newfile
}

################################################################################
