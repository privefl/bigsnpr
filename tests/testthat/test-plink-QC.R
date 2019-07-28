################################################################################

context("PLINK_QC")

################################################################################

# Get PLINK executable
PLINK <- download_plink()

unlink(plink2 <- download_plink2(AVX2 = FALSE))
expect_false(grepl(plink2, "avx2"))

PLINK2 <- download_plink2()

################################################################################

# Filter on MAF
library(magrittr)
bedfile <- snp_attachExtdata() %>%
  snp_writeBed(tempfile(fileext = ".bed"))

snp_plinkQC(plink.path = PLINK,
            prefix.in = sub_bed(bedfile),
            prefix.out = tempfile(),
            maf = 0.2,
            extra.options = "--allow-no-sex") %>%
  snp_readBed() %>%
  snp_attach() %>%
  extract2("genotypes") %>%
  snp_MAF() %>%
  is_greater_than(0.2) %>%
  all() %>%
  expect_true()

################################################################################

# IBD
df.pair <- snp_plinkIBDQC(plink.path = PLINK,
                          bedfile.in = bedfile,
                          bedfile.out = tempfile(fileext = ".bed"),
                          pi.hat = 0.1,
                          do.blind.QC = FALSE,
                          extra.options = "--allow-no-sex")
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
  snp_plinkRmSamples(plink.path = PLINK,
                     bedfile.in = bedfile,
                     bedfile.out = tempfile(fileext = ".bed"),
                     df.or.files = df.pair) %>%
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

bedfile2 <- snp_plinkKINGQC(PLINK2, bedfile)
expect_identical(readLines(sub_bed(bedfile2, ".king.cutoff.out.id")),
                 c("#FID\tIID", "POP1\tIND2"))

bedfile3 <- snp_plinkKINGQC(PLINK2, bedfile, thr.king = 0.3,
                            bedfile.out = tempfile(fileext = ".bed"))
expect_identical(readLines(sub_bed(bedfile3, ".king.cutoff.out.id")), "#FID\tIID")

expect_error(snp_plinkKINGQC(PLINK, bedfile), "This requires PLINK v2")

################################################################################
