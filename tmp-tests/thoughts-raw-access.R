x <- sample.int(1e6)

x2 <- x - 1
x3 <- x2 %% 4
x4 <- x2 %/% 4

print(all.equal(x4 * 4 + x3 + 1, x))

x <- integer(10)
x[] <- rep(1:3, 10)


# match each possible code
getCode <- function() {
  all.raws <- as.raw(0:255)
  geno.raw <- as.logical(rawToBits(all.raws))
  s <- c(TRUE, FALSE)
  geno1 <- geno.raw[s]
  geno2 <- geno.raw[!s]
  geno <- geno1 + geno2
  geno[geno1 & !geno2] <- -128
  dim(geno) <- c(4, 256)
  geno
}
geno <- getCode()

`[.bigBed` <- function(x) {
  x$X
}
