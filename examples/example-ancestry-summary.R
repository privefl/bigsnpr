## Toy example

X <- matrix(rbeta(4000, 1, 3), ncol = 4)
y <- X %*% c(0.2, 0.1, 0.6, 0.1)
res <- snp_ancestry_summary(y, X)
str(res)


## Real GWAS summary statistics: Epilepsy (supposedly in EUR+EAS+AFR)

\dontrun{

library(dplyr)

gz <- runonce::download_file(
  "http://www.epigad.org/gwas_ilae2018_16loci/all_epilepsy_METAL.gz",
  dir = "tmp-data")
readLines(gz, n = 3)
sumstats <- bigreadr::fread2(
  gz, select = c("CHR", "BP", "Allele2", "Allele1", "Freq1"),
  col.names = c("chr", "pos", "a0", "a1", "freq")
) %>%
  mutate_at(3:4, toupper)

all_freq <- bigreadr::fread2(
  runonce::download_file("https://figshare.com/ndownloader/files/31218445",
                         dir = "tmp-data", fname = "ref_freqs.csv.gz"))

matched <- bigsnpr::snp_match(
  mutate(sumstats, chr = as.integer(chr), beta = 1),
  all_freq[1:5], return_flip_and_rev = TRUE
) %>%
  mutate(freq = ifelse(`_REV_`, 1 - freq, freq))

X <- as.matrix(all_freq[matched$`_NUM_ID_`, -(1:5)])
y <- matched$freq
snp_ancestry_summary(y, X)
}
