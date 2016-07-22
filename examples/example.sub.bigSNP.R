\dontrun{

  bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")

  # Creating directory for backing files
  if (!dir.exists("backingfiles")) dir.create("backingfiles")

  if (!file.exists("backingfiles/test_doc.bk")) {
    # Reading the bedfile and storing the data
    test <- BedToBig(bedfile, 50, "test_doc")
  } else {
    # Loading it from backing files
    test <- AttachBigSNP("test_doc")
  }

  # dimensions
  print(rbind(dim(test$genotypes), dim(test$fam), dim(test$map)))

  # taking only the first 50 individuals and 500 SNPs at random
  ind.row <- 1:50
  ind.col <- sort(sample(ncol(test$genotypes), 500))
  test2 <- sub.bigSNP(test, ind.row, ind.col)

  # dimensions
  print(rbind(dim(test2$genotypes), dim(test2$fam), dim(test2$map)))
  # new backing files
  print(test2$backingfile)
  print(test2$backingpath)


  # removing the 100th and 120th individuals
  ind.row.del <- c(100, 120)
  test3 <- sub.bigSNP(test, -ind.row.del)

  # dimensions
  print(rbind(dim(test3$genotypes), dim(test3$fam), dim(test3$map)))
  print(test3$fam$sample.ID[98:120])
  # new backing files
  print(test3$backingfile) # incrementation of the file's number
}
