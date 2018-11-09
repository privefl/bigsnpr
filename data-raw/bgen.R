bgen_file <- "inst/testdata/example.8bits.bgen"
system(glue::glue("bgenix -g {bgen_file} -index"))

library(dplyr)
db_con <- RSQLite::dbConnect(RSQLite::SQLite(), paste0(bgen_file, ".bgi"))
infos <- tbl(db_con, "Variant") %>%
 collect()
RSQLite::dbDisconnect(db_con)

res <- rbgen::bgen.load(filename = bgen_file, rsids = infos$rsid)
str(res)
stopifnot(all(res$ploidy == 2))
stopifnot(any(res$phased) == FALSE)

variants <- mutate_if(res$variants, is.factor, as.character)
saveRDS(as_tibble(setNames(variants[c(1, 3, 2, 5:6)], bigsnpr:::NAMES.MAP[-3])),
        "inst/testdata/bgen_variants.rds")

# compute dosages from probabilities
DIM <- dim(res$data)
dim(res$data) <- c(prod(DIM[1:2]), DIM[3])
dosages <- res$data %*% 0:2
dim(dosages) <- DIM[1:2]
saveRDS(t(dosages), "inst/testdata/bgen_dosages.rds")
