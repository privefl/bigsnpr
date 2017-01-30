bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
print(bedfile)

# Reading the bedfile and storing the data in directory "backingfiles"
PATH <- "backingfiles/test_doc.bk"
if (!file.exists(PATH))
  PATH <- snp_readBed(bedfile, backingfile = "test_doc")

# Loading the data from backing files
test <- snp_attach(PATH)

str(test)
print(dim(test$genotypes))
print(test$genotypes[1:8, 1:8])
