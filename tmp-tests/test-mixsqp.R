library(dplyr)

all_freq <- bigreadr::fread2(
  runonce::download_file("https://figshare.com/ndownloader/files/31218445",
                         dir = "tmp-data", fname = "ref_freqs.csv.gz"))

# all_freq[sample(nrow(all_freq), 2), -(1:5)]


# Epilepsy
gz <- runonce::download_file(
  "http://www.epigad.org/gwas_ilae2018_16loci/all_epilepsy_METAL.gz",
  dir = "tmp-data")
readLines(gz, n = 3)
sumstats <- bigreadr::fread2(
  gz, select = c("CHR", "BP", "Allele2", "Allele1", "Freq1"),
  col.names = c("chr", "pos", "a0", "a1", "freq")
) %>%
  mutate_at(3:4, toupper)

# PAGE
gz <- runonce::download_file(
  "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST008001-GCST009000/GCST008053/WojcikG_PMID_invn_rheight_alls.gz",
  dir = "tmp-data")
sumstats <- bigreadr::fread2(
  gz, select = c("Chr", "Position_hg19", "Other-allele", "Effect-allele", "Effect-allele-frequency"),
  col.names = c("chr", "pos", "a0", "a1", "freq"))

pop <- c("African-American" = 17286, "Hispanic-Latino" = 22192, "Asian" = 4680,
         "Native Hawaiian" = 3939, "Native American" = 647,	Other	= 1052)
round(sort(pop / sum(pop), decreasing = TRUE), 3)
# Hispanic-Latino African-American    Asian  Native Hawaiian    Other  Native American
#           0.446            0.347    0.094            0.079    0.021            0.013


matched <- bigsnpr::snp_match(
  mutate(sumstats, chr = as.integer(chr), beta = 1),
  all_freq[1:5], return_flip_and_rev = TRUE
) %>%
  mutate(freq = ifelse(`_REV_`, 1 - freq, freq))

X <- as.matrix(all_freq[matched$`_NUM_ID_`, -(1:5)])
y <- matched$freq


res0 <- bigsnpr::snp_ancestry_summary(y, X)

all_N <- bigreadr::fread2("../freq-ancestry/all_prop.csv")$N

logL <- do.call("cbind", lapply(seq_along(all_N), function(k) {
  N <- all_N[k]
  p <- X[, k]
  dbeta(y, 1 + 2 * N * p, 1 + 2 * N * (1 - p), log = TRUE)
}))
res <- mixsqp::mixsqp(logL, log = TRUE)

round(res$x, 4)

round(rbind(res$x, res0$coef), 4)
