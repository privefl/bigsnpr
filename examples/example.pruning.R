# Creating directory for backing files
if (!dir.exists("backingfiles")) dir.create("backingfiles")

if (!file.exists("backingfiles/test_doc.bk")) {
  # Reading the bedfile and storing the data
  bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
  test <- BedToBig(bedfile, 50, "test_doc")
} else {
  # Loading it from backing files
  test <- AttachBigSNP("test_doc")
}

# generating random phenotypes
test$fam$pheno <- sample(c(-1, 1), size = nrow(test$fam),
                         replace = TRUE)

ind <- Prune(test)
# ind <- Prune(test, ncores = 2)

R2 <- RsqClass(test$genotypes, test$fam$pheno)

par.save <- par(mfrow = c(1, 2))

plot(R2, type = 'h')
plot(ind, R2[ind], type = "h")

par(par.save)
