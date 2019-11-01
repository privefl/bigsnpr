################################################################################

context("MATCH")

################################################################################

sumstats <- data.frame(
  chr = 1,
  pos = c(86303, 86331, 162463, 752566, 755890, 758144),
  a0 = c("T", "G", "C", "A", "T", "G"),
  a1 = c("G", "A", "T", "G", "A", "A"),
  beta = c(-1.868, 0.250, -0.671, 2.112, 0.239, 1.272),
  p = c(0.860, 0.346, 0.900, 0.456, 0.776, 0.383)
)

info_snp <- data.frame(
  chr = 1,
  rsid = c("rs2949417", "rs115209712", "rs143399298", "rs3094315", "rs3115858"),
  a0 = c("T", "A", "G", "A", "T"),
  a1 = c("G", "G", "A", "G", "A"),
  pos = c(86303, 86331, 162463, 752566, 755890)
)

expect_message(matched1 <- snp_match(sumstats, info_snp),
               "4 variants have been matched; 1 were flipped and 1 were reversed.")
expect_equal(dim(matched1), c(4, 9))
expect_equal(matched1$beta, sumstats$beta[1:4] * c(1, -1, 1, 1))
expect_message(matched2 <- snp_match(sumstats, info_snp, strand_flip = FALSE),
               "4 variants have been matched; 0 were flipped and 1 were reversed.")
expect_equal(dim(matched2), c(4, 9))
expect_equal(matched2$beta, sumstats$beta[c(1:2, 4:5)] * c(1, -1, 1, 1))


sumstats2 <- data.frame(
  chr = 1,
  rsid = c("rs2949417", "rs115209712", "rs143399298", "rs3094315", "rs3115858", NA),
  a0 = c("T", "G", "C", "A", "T", "G"),
  a1 = c("G", "A", "T", "G", "A", "A"),
  pos = 10 + c(86303, 86331, 162463, 752566, 755890, 758144),
  beta = 1
)

expect_error(snp_match(sumstats2[-6], info_snp), "'chr, pos, a0, a1, beta'")
expect_error(snp_match(sumstats2, info_snp[-5]), "'chr, pos, a0, a1'")
expect_error(snp_match(sumstats2, info_snp),
             "Not enough variants have been matched.")
expect_equal(dim(snp_match(sumstats2, info_snp, match.min.prop = 0)), c(0, 9))
snp_info <- snp_match(sumstats2, info_snp, join_by_pos = FALSE)
expect_equal(dim(snp_info), c(4, 9))
expect_equal(snp_info$beta, c(1, -1, 1, 1))
expect_equal(snp_info$pos.ss, snp_info$pos + 10)
expect_equal(snp_info[c(1:4, 8)], info_snp[snp_info$`_NUM_ID_`, ],
             check.attributes = FALSE)
expect_equal(snp_info[c(1, 2, 5)], sumstats2[snp_info$`_NUM_ID_.ss`, c(1, 2, 5)],
             check.attributes = FALSE)

################################################################################

expect_equal(same_ref(ref1 = info_snp$a1, alt1 = info_snp$a0,
                      ref2 = sumstats$a1[1:5], alt2 = sumstats$a0[1:5]),
             c(TRUE, FALSE, TRUE, TRUE, TRUE))

################################################################################
