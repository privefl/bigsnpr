setwd("../bigsnpr")

library(dplyr)

all_freq <- bigreadr::fread2(
  runonce::download_file("https://figshare.com/ndownloader/files/31218445",
                         dir = "tmp-data", fname = "ref_freqs.csv.gz"))

# T2D from FinnGen
(N <- 29193 + 182573)
sumstats <- bigreadr::fread2("../freq-ancestry//tmp-data/summary_stats_finngen_R5_T2D.gz") %>%
  transmute(rsid = rsids, chr = `#chrom`, pos, a0 = ref, a1 = alt,
            freq = (2 * (n_hom_cases + n_hom_controls) +
                      (n_het_cases + n_het_controls)) / (2 * N)) %>%
  filter(chr %in% 1:22)


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
  all_freq[1:5], return_flip_and_rev = TRUE, join_by_pos = FALSE
) %>%
  mutate(freq = ifelse(`_REV_`, 1 - freq, freq))

X <- as.matrix(all_freq[matched$`_NUM_ID_`, -(1:5)])
y <- matched$freq

D <- crossprod(X)
eig_min <- min(eigen(D, symmetric = TRUE, only.values = TRUE)$values)

get_res <- function(coef) {
  diag(D) <- diag(D) - coef * eig_min
  res <- quadprog::solve.QP(
    Dmat = D,
    dvec = crossprod(y, X),
    Amat = cbind(1, diag(ncol(X))),
    bvec = c(1, rep(0, ncol(X))),
    meq  = 1
  )$solution
  round(100 * setNames(res, colnames(X)), 1)
}

(SEQ <- setNames(nm = c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99)))

sapply(SEQ, get_res) %>%
  ifelse(. == 0, "", .) %>%
  xtable::xtable(align = "|l|c|c|c|c|c|c|c|c|") %>%
  print(caption.placement = "top", hline.after = c(-1, 0, 4, 5, 10, 11, 14, 17),
        include.rownames = TRUE, file = "../freq-ancestry/ancestry-proportions.tex")
