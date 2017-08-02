################################################################################

context("PLINK_QC")

################################################################################

# Get PLINK executable
download.file("https://www.cog-genomics.org/static/bin/plink170725/plink_linux_x86_64.zip",
              destfile = (plink.zip <- tempfile()))
PLINK <- unzip(plink.zip, files = "plink", exdir = tempdir())
Sys.chmod(PLINK, mode = (file.info(PLINK)$mode | "111"))

################################################################################

# Filter on MAF
library(magrittr)
bedfile <- snp_attachExtdata() %>%
  snp_writeBed(tempfile(fileext = ".bed"))

snp_plinkQC(plink.path = PLINK,
            prefix.in = sub("\\.bed$", "", bedfile),
            prefix.out = tempfile(),
            maf = 0.2,
            extra.options = "--allow-no-sex") %>%
  snp_attachExtdata() %>%
  extract2("genotypes") %>%
  snp_MAF() %>%
  magrittr::is_greater_than(0.2) %>%
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
  attach.BM() %>%
  as.matrix() %>%
  t() %>%
  cor()

expect_lt(t.test(K[ind.pair], K, alternative = "greater")$p.value, 1e-16)

indiv.keep <-
  snp_plinkRmSamples(plink.path = PLINK,
                     bedfile.in = bedfile,
                     bedfile.out = tempfile(fileext = ".bed"),
                     df.or.files = df.pair) %>%
  snp_attachExtdata() %>%
  extract2("fam") %>%
  extract2("sample.ID")

expect_length(union(indiv.keep, df.pair$IID1), nrow(K))
expect_length(intersect(indiv.keep, df.pair$IID1), 0)

################################################################################
