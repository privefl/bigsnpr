################################################################################

context("PLINK_QC")

################################################################################

# Get PLINK executable
plink <- download_plink()
expect_false(grepl("plink2", plink))

# plink2 <- download_plink2(overwrite = TRUE)
# expect_true(grepl("plink2", plink2))
# plink2.version <- system(paste(plink2, "--version"), intern = TRUE)
# expect_true(grepl("AVX2", plink2.version))

plink2 <- download_plink2(overwrite = TRUE, AVX2 = FALSE)
expect_true(grepl("plink2", plink2))
plink2.version <- system(paste(plink2, "--version"), intern = TRUE)
expect_false(grepl("AVX2", plink2.version))

################################################################################

# Filter on MAF
library(magrittr)
bedfile <- snp_attachExtdata() %>%
  snp_writeBed(tempfile(fileext = ".bed"))

snp_plinkQC(plink.path = plink,
            prefix.in = sub_bed(bedfile),
            prefix.out = tempfile(),
            maf = 0.2,
            extra.options = "--allow-no-sex",
            verbose = FALSE) %>%
  snp_readBed() %>%
  snp_attach() %>%
  extract2("genotypes") %>%
  snp_MAF() %>%
  is_greater_than(0.2) %>%
  all() %>%
  expect_true()

################################################################################

# IBD
df.pair <- snp_plinkIBDQC(plink.path = plink,
                          bedfile.in = bedfile,
                          bedfile.out = tempfile(fileext = ".bed"),
                          pi.hat = 0.1,
                          do.blind.QC = FALSE,
                          extra.options = "--allow-no-sex",
                          verbose = FALSE)
ind.pair <- df.pair %>%
  extract(c("IID1", "IID2")) %>%
  sapply(function(x) as.numeric(gsub("IND", "", x)) + 1)

K <- snp_attachExtdata() %>%
  extract2("genotypes") %>%
  extract(,) %>%
  t() %>%
  cor()

expect_lt(t.test(K[ind.pair], K, alternative = "greater")$p.value, 1e-16)

indiv.keep <-
  snp_plinkRmSamples(plink.path = plink,
                     bedfile.in = bedfile,
                     bedfile.out = tempfile(fileext = ".bed"),
                     df.or.files = df.pair,
                     verbose = FALSE) %>%
  snp_readBed() %>%
  snp_attach() %>%
  extract2("fam") %>%
  extract2("sample.ID")

expect_length(union(indiv.keep, df.pair$IID1), nrow(K))
expect_length(intersect(indiv.keep, df.pair$IID1), 0)

################################################################################

obj.snp <- snp_attachExtdata()
G <- obj.snp$genotypes
G[1, ] <- round(colMeans(G[2:3, ]))
bedfile <- snp_writeBed(obj.snp, tempfile(fileext = ".bed"))

expect_error(snp_plinkKINGQC(plink, bedfile), "This requires PLINK v2")

if (.Machine$sizeof.pointer == 8) {
  bedfile2 <- snp_plinkKINGQC(plink2, bedfile, thr.king = 0.177, verbose = FALSE)
  expect_identical(readLines(sub_bed(bedfile2, ".king.cutoff.out.id")),
                   c("#FID\tIID", "POP1\tIND2"))
  keep <- snp_plinkKINGQC(plink2, bedfile, thr.king = 0.177,
                          make.bed = FALSE, verbose = FALSE)
  expect_false(keep[3])
  expect_true(all(keep[-3]))

  bedfile3 <- snp_plinkKINGQC(plink2, bedfile, thr.king = 0.3,
                              bedfile.out = tempfile(fileext = ".bed"),
                              verbose = FALSE)
  expect_identical(readLines(sub_bed(bedfile3, ".king.cutoff.out.id")), "#FID\tIID")
  keep2 <- snp_plinkKINGQC(plink2, bedfile, thr.king = 0.3,
                           make.bed = FALSE, verbose = FALSE)
  expect_true(all(keep2))
}

unlink(plink2)

################################################################################
