(bedfile <- system.file("extdata", "example.bed", package = "bigsnpr"))

# Reading the bedfile and storing the data in directory "backingfiles"
rds <- snp_readBed(bedfile, backingfile = "test_doc")

# Loading the data from backing files
test <- snp_attach(rds)

str(test)
dim(G <- test$genotypes)
attach.BM(G)[1:8, 1:8]



# cleaning
bk <- sub("\\.rds$", ".bk", rds)
file.remove(rds, bk)
