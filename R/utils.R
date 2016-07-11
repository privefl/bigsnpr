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

LimsChr <- function(infos, qc = FALSE) {
  if (qc) {
    map <- infos$map$chromosome[-infos$indQC]
  } else {
    map <- infos$map$chromosome
  }

  upper <- cumsum(rle(map)$length)
  lower <- c(1, upper[-length(upper)] + 1)
  cbind(lower, upper)
}

################################################################################
