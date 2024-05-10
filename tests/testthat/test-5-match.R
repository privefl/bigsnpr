################################################################################

context("MATCH")

# test_that()

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
expect_message(snp_match(data.table::as.data.table(sumstats), info_snp),
               "4 variants have been matched; 1 were flipped and 1 were reversed.")
expect_message(snp_match(sumstats, data.table::as.data.table(info_snp)),
               "4 variants have been matched; 1 were flipped and 1 were reversed.")

expect_message(matched3 <- snp_match(sumstats, info_snp, return_flip_and_rev = TRUE),
               "4 variants have been matched; 1 were flipped and 1 were reversed.")
expect_equal(dim(matched3), c(4, 11))
expect_equal(matched3[["_FLIP_"]], c(FALSE, FALSE,  TRUE, FALSE))
expect_equal(matched3[["_REV_"]],  c(FALSE,  TRUE, FALSE, FALSE))


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
expect_error(snp_match(sumstats2, info_snp), "No variant has been matched.")
snp_info <- snp_match(sumstats2, info_snp, join_by_pos = FALSE)
expect_equal(dim(snp_info), c(4, 9))
expect_equal(snp_info$beta, c(1, -1, 1, 1))
expect_equal(snp_info$pos.ss, snp_info$pos + 10)
expect_equal(snp_info[c(1:4, 8)], info_snp[snp_info$`_NUM_ID_`, ],
             check.attributes = FALSE)
expect_equal(snp_info[c(1, 2, 5)], sumstats2[snp_info$`_NUM_ID_.ss`, c(1, 2, 5)],
             check.attributes = FALSE)
expect_message(
  snp_info <- snp_match(sumstats2[c(1, 1:6), ], info_snp, join_by_pos = FALSE),
  "Some duplicates were removed.")
expect_equal(dim(snp_info), c(3, 9))

################################################################################

expect_equal(same_ref(ref1 = info_snp$a1, alt1 = info_snp$a0,
                      ref2 = sumstats$a1[1:5], alt2 = sumstats$a0[1:5]),
             c(TRUE, FALSE, TRUE, TRUE, TRUE))

################################################################################

test_that("snp_asGeneticPos() works", {

  info <- data.frame(
    chr = rep(1L, 5L),
    pos = c(853954L, 854250L, 864938L, 870645L, 873558L),
    rsid = c("rs1806509", "rs7537756", "rs2340587", "rs28576697", "rs1110052"))

  map <- data.frame(
    V1 = c("rs1806509", "rs7537756", "rs28576697", "rs1110052"),
    V2 = c(853954L, 854250L, 870645L, 873558L),
    V3 = c(0.194323402834, 0.194576977815, 0.202835640491, 0.203874368612))
  bigreadr::fwrite2(map, file.path(tempdir(), "chr1.OMNI.interpolated_genetic_map"))

  res1 <- snp_asGeneticPos(info$chr, info$pos, dir = tempdir(), ncores = 2)
  expect_equal(res1, map$V3[c(1, 2, 3, 3, 4)])

  res2 <- snp_asGeneticPos(info$chr, info$pos, dir = tempdir(), ncores = 2,
                           rsid = info$rsid)
  expect_equal(res2[-3], map$V3)
  expect_gt(res2[3], res2[2])
  expect_lt(res2[3], res2[4])

  skip_if(is_cran)
  skip_on_covr()
  skip_if_offline("raw.githubusercontent.com")

  info2 <- data.frame(chr = 22, pos = 18206376)
  res3 <- snp_asGeneticPos(info2$chr, info2$pos, dir = tempdir())
  expect_equal(res3, 3.682628, tolerance = 1e-5)
  res4 <- snp_asGeneticPos(info2$chr, info2$pos, dir = tempdir(), type = "hapmap")
  expect_equal(res4, 5.905713, tolerance = 1e-5)
})

################################################################################

test_that("snp_asGeneticPos2() works", {

  info <- data.frame(
    chr = rep(1L, 5L),
    pos = c(853954L, 854250L, 864938L, 870645L, 873558L))

  map <- data.frame(
    chr = 1,
    pos = c(853954L, 854250L, 870645L, 873558L),
    pos_cM = c(0.194323402834, 0.194576977815, 0.202835640491, 0.203874368612))

  res1 <- snp_asGeneticPos2(info$chr, info$pos, genetic_map = map)
  expect_equal(res1[-3], map$pos_cM[c(1, 2, 3, 4)])
  expect_equal(res1[3],  map$pos_cM[3], tolerance = 0.01)

  # if outside limits, just use boundaries
  res2 <- snp_asGeneticPos2(c(1, 1), c(0, 1e6), genetic_map = map)
  expect_equal(res2, range(map$pos_cM))

  # some errors
  expect_error(snp_asGeneticPos2(c(1, 1), c(0, 1e6), genetic_map = map[1:2]),
               "'genetic_map' should have element 'pos_cM'.")
  expect_error(snp_asGeneticPos2(c(2, 2), c(0, 1e6), genetic_map = map),
               "Chromosome '2' not found in `genetic_map`.")

})

################################################################################

test_that("snp_ancestry_summary() works (with no projection here)", {

  X <- matrix(rbeta(4000, 1, 3), ncol = 4)
  prop <- c(0.2, 0.1, 0.6, 0.1)
  y <- X %*% prop
  res <- snp_ancestry_summary(y, X, Matrix::Diagonal(nrow(X), 1), rep(1, nrow(X)))
  expect_equal(res, prop, check.attributes = FALSE)
  expect_equal(attr(res, "cor_pred"), 1)

  expect_warning(res2 <- snp_ancestry_summary(
    y, X[, -1], Matrix::Diagonal(nrow(X), 1), rep(1, nrow(X))),
    "The solution does not perfectly match the frequencies.")
  expect_true(all(res2 > prop[-1]))
  expect_equal(res2, prop[-1], tolerance = 0.1, check.attributes = FALSE)
  expect_equal(attr(res2, "cor_pred"), drop(cor(y, X[, -1] %*% res2)))

  expect_error(res3 <- snp_ancestry_summary(
    y, X[, -3], Matrix::Diagonal(nrow(X), 1), rep(1, nrow(X)), min_cor = 0.6),
    "Correlation between frequencies is too low:")
})

################################################################################
