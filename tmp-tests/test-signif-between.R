x <- rnorm(1e7)
x2 <- sort(x, decreasing = TRUE)
a <- x2[6e6]
b <- x2[6e6+1]

signif(a)
signif(b)

i <- 1
repeat {
  s <- signif(b, digits = i)
  if ((s >= b) && (a > s)) {
    break
  } else {
    if (i == 8) break
    i <- i + 1
  }
}


