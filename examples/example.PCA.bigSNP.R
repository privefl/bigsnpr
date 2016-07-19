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

  # constructing a bigSNP with 517 rows and 227,100 columns
  # Ok, "sub" is not well adapted when using it like this..
  test2 <- sub.bigSNP(test, ind.col = rep(1:ncol(test$genotypes), times = 50))
  print(dim(test2$genotypes))

  # Computing with Eigen library: 8 seconds on my computer
  print(system.time(test3 <- BigXYt(test2, 1000)))
  print(test3[1:10, 1:5])

  # Computing with base R: 19 seconds on my computer
  # Less than 5 seconds when using Microsoft R Open
  print(system.time(test4 <- BigXYt(test2, 1000, use.Eigen = FALSE)))
  print(test4[1:10, 1:5])
}
