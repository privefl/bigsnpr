printf <- function(...) cat(sprintf(...))


CutBySize <- function(m, block.size) {
  nb = ceiling(m / block.size)
  int = m / nb

  upper = round(1:nb * int)
  lower = c(1, upper[-nb] + 1)
  size = c(upper[1], diff(upper))

  return( cbind(lower, upper, size) )
}
