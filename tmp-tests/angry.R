a <- big.matrix(4, 4, type = "char")

a[] <- sample(c(0, 1, 2, NA), size = 16, replace = TRUE)

a[]

m <- big.matrix(10, 1, type = 'short')
m[4,1] <- as.raw(130)
