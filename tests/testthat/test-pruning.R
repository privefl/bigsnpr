################################################################################

context("PRUNING")

test <- snp_attachExtdata()
G <- test$genotypes
ids <- test$map$marker.ID
chrs <- test$map$chromosome
pos <- test$map$physical.pos

expect_length(snp_indLRLDR(chrs, pos), 0)
pos2 <- round(runif(length(pos), 48060567, 52060567))
expect_error(snp_indLRLDR(chrs, pos2))
expect_length(snp_indLRLDR(chrs, sort(pos2)), ncol(G))

################################################################################

expect_warning(ind.keep <- snp_pruning(G, chrs), "Pruning is deprecated")

ind.true <- scan(
  file = system.file("testdata", "pruning_ct.prune.in", package = "bigsnpr"),
  what = "",
  quiet = TRUE
)

test_that("Same indices when pruning", {
  expect_equal(ids[ind.keep], ind.true)
})

################################################################################

expect_warning(
  ind.keep <- snp_pruning(G, chrs, size = 200, infos.pos = pos, thr.r2 = 0.1),
  "Pruning is deprecated")

ind.true <- scan(
  file = system.file("testdata", "pruning_kb.prune.in", package = "bigsnpr"),
  what = "",
  quiet = TRUE
)

test_that("Same indices when pruning with size in kb", {
  expect_gte(length(ind.keep), length(ind.true))
})

################################################################################
