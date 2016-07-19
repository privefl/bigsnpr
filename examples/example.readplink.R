\dontrun{

  bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")

  # Creating directory for backing files
  if (!dir.exists("backingfiles")) dir.create("backingfiles")

  # Reading the bedfile and storing the data in directory "backingfiles"
  if (!file.exists("backingfiles/test_doc.bk"))
    test <- BedToBig(bedfile, 50, "test_doc")

  # Removing the R object
  rm(test)

  # Loading it from backing files
  test <- AttachBigSNP("test_doc")

  str(test)
  print(dim(test$genotypes))
  print(test$genotypes[1:8, 1:8])
}
