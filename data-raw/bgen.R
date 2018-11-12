# https://bitbucket.org/gavinband/bgen/src/default/example/example.8bits.bgen
bgen_file <- "tmp-data/example.8bits.bgen"
# system(glue::glue("bgenix -g {bgen_file} -index"))

saveRDS(readBin(bgen_file, what = raw(), n = 1e6), "inst/testdata/bgen_example.rds")
saveRDS(readBin(paste0(bgen_file, ".bgi"), what = raw(), n = 1e6),
        "inst/testdata/bgi_example.rds")

library(dplyr)
db_con <- RSQLite::dbConnect(RSQLite::SQLite(), paste0(bgen_file, ".bgi"))
infos <- tbl(db_con, "Variant") %>%
 collect()
RSQLite::dbDisconnect(db_con)

res <- rbgen::bgen.load(filename = bgen_file, rsids = infos$rsid)
str(res)
stopifnot(all(res$ploidy == 2))
stopifnot(any(res$phased) == FALSE)

variants <- res$variants %>%
  as_tibble() %>%
  mutate_if(is.factor, as.character) %>%
  mutate(marker.ID = sub("^RSID_", "SNPID_", rsid), allele2 = allele1) %>%
  select(chromosome, marker.ID, rsid, physical.pos = position,
         allele1 = allele0, allele2)

saveRDS(variants, "inst/testdata/bgen_variants.rds")

# compute dosages from probabilities
DIM <- dim(res$data)
dim(res$data) <- c(prod(DIM[1:2]), DIM[3])
dosages <- res$data %*% 0:2
dim(dosages) <- DIM[1:2]
saveRDS(t(dosages), "inst/testdata/bgen_dosages.rds")
