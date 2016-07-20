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

  # create some ramdon class
  test$fam$affection <- sample(c(-1, 1),
                               length(test$fam$affection),
                               replace = TRUE)

  # Get predictions from half the data
  n <- nrow(test$genotypes)
  ind.train <- sort(sample(n, n/2))
  test2 <- GBLUP(test, block.size = 1000, ind.train = ind.train)
  print(AucSampleConf(test2, test$fam$affection[-ind.train], nboot = 1e4))
}
