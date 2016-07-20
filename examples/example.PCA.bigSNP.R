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

  # Get the first 10 PCs
  test2 <- PCA.bigSNP(test, block.size = 1000, k = 10)

  # Make a nice plot out of the two first PCs
  grp.pop <- test$fam$family.ID
  grp.pop[grp.pop %in% paste0("POP", 1:5)] <- 1
  grp.pop[grp.pop %in% paste0("POP", 6:11)] <- 2
  grp.pop[grp.pop %in% paste0("POP", 12:19)] <- 3
  x <- test2[, 1]
  is.right = (max(x) > max(-x))
  y <- test2[, 2]
  is.top = (max(y) > max(-y))
  legend.pos <- ifelse(is.top, ifelse(is.right, "topright", "topleft"),
                       ifelse(is.right, "bottomright", "bottomleft"))
  plot(x, y, xlab = "PC1", ylab = "PC2", col = grp.pop)
  legend(legend.pos, legend = c("POP1 to POP5",
                                "POP6 to POP11",
                                "POP12 to POP19"),
         col = 1:3, pch = 1, cex = 0.8)

  # Train PCs with only half the data
  n <- nrow(test$genotypes)
  ind.train <- sort(sample(n, n/2))
  test3 <- PCA.bigSNP(test, block.size = 1000, k = 10,
                      ind.train = ind.train)

  # Compare with prcomp
  test4 <- prcomp(test$genotypes[,], scale. = TRUE)
  plot(test4$x, col = grp.pop)
}
