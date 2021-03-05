library(magrittr)
bgen_file <- "~/test.bgen"
system.file("testdata", "bgen_example.rds", package = "bigsnpr") %>%
  readRDS() %>% writeBin(bgen_file, useBytes = TRUE)
system.file("testdata", "bgi_example.rds",  package = "bigsnpr") %>%
  readRDS() %>% writeBin(paste0(bgen_file, ".bgi"), useBytes = TRUE)

variants <- readRDS(system.file("testdata", "bgen_variants.rds", package = "bigsnpr"))
IDs <- with(variants, paste(1, physical.pos, allele1, allele2, sep = "_"))

test <- bigsnpr::snp_readBGEN(bgen_file, tempfile(), list(IDs))

