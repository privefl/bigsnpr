source('R/utils.R')
Rcpp::sourceCpp('tmp-tests/fstream.cpp')

tab <- getInverseCode()

require(bigsnpr)
celiac <- AttachBigSNP("celiac")
X <- celiac$genotypes

test <- writebina("test.bed", X@address, tab)

bedfile <- "test.bed"
bed <- file(bedfile, open = "rb")
test <- readBin(bed, "raw", 1e6)
close(bed)

bedfile2 <- "../Dubois2010_data/FinnuncorrNLITUK1UK3hap300.bed"
bed2 <- file(bedfile2, open = "rb")
test2 <- readBin(bed2, "raw", 1e6)
close(bed2)
all.equal(test, test2) # OK
