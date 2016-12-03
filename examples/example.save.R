# Creating directory for backing files
if (!dir.exists("backingfiles")) dir.create("backingfiles")

if (!file.exists("backingfiles/test_doc.bk")) {
  # Reading the bedfile and storing the data
  bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
  test <- BedToBig(bedfile, backingpath = "test_doc")
} else {
  # Loading it from backing files
  test <- AttachBigSNP("test_doc")
}

# Just after reading
print(rle(test$fam$family.ID))
print(test$fam$pop) # NULL

# Get populations from slot "family.ID"
test <- GetPops(test)
print(rle(test$fam$pop))

# Get populations clusters from external files
files <- system.file("extdata", paste0("cluster", 1:3), package = "bigsnpr")
data.table::fread(files[1])
test <- GetPops(x = test, pop.files = files, col.sample.ID = 2, col.family.ID = 3)
print(rle(test$fam$pop))
