
library(bigsnpr)

bedfile <- "../Dubois2010_data/FinnuncorrNLITUK1UK3hap300.bed"

system.time(
  test1 <- snp_attach(snp_readBed(bedfile, tempfile()))
)

system.time({
  backingfile <- tempfile()

  # Get other files
  bimfile <- sub("\\.bed$", ".bim", bedfile)
  famfile <- sub("\\.bed$", ".fam", bedfile)

  # Read map and family files
  fam <- data.table::fread(famfile, data.table = FALSE)
  # names(fam) <- NAMES.FAM
  bim <- data.table::fread(bimfile, data.table = FALSE)
  # names(bim) <- NAMES.MAP

  n <- nrow(fam)
  m <- nrow(bim)

  # Prepare Filebacked Big Matrix
  bigGeno <- FBM.code256(
    nrow = n,
    ncol = m,
    code = bigsnpr:::CODE_012,
    backingfile = backingfile,
    init = NULL,
    create_bk = TRUE,
    save = FALSE
  )

  # Fill the FBM from bedfile
  # Rcpp::sourceCpp('tmp-tests/test-read-bed-mmap.cpp')
  fillMat(pcadapt:::bedXPtr(bedfile, n, m), bigGeno)

  # Create the bigSNP object
  snp.list <- structure(list(genotypes = bigGeno,
                             fam = fam,
                             map = bim),
                        class = "bigSNP")

  # save it and return the path of the saved object
  rds <- sub("\\.bk$", ".rds", bigGeno$backingfile)
  saveRDS(snp.list, rds)
  test2 <- snp_attach(rds)
})


test1$genotypes[1:5, 1:5]
test2$genotypes[1:5, 1:5]


# snp_writeBed(test1, bedfile = "test.bed", ind.row = rep(seq_len(n), 10))
