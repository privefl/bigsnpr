require(bigsnpr)

celiac <- AttachBigSNP("celiac")
X <- celiac$genotypes
n <- nrow(X)
complete <- -n %% 4
bsz <- ceiling(n/4)
x <- X[, 1:10]
m2 <- ncol(x)
x2 <- rbind(replace(x, is.na(x), 3L), matrix(0L, complete, m2))

r <- getInverseCode()
dim(x2) <- c(4, bsz * m2)
x3 <- t(x2 + 1L)
test <- r[x3]
head(test)

bedfile <- "../Dubois2010_data/FinnuncorrNLITUK1UK3hap300.bed"
bed <- file(bedfile, open = "rb")
magic <- readBin(bed, "raw", 3)
test2 <- readBin(bed, "raw", length(test))
close(bed)
all.equal(test, test2) # OK

