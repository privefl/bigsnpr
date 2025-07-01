################################################################################

context("FST")

skip_on_os("solaris")
skip_if(is_cran)
skip_if_offline("www.cog-genomics.org")
skip_if_offline("s3.amazonaws.com")

plink <- download_plink(verbose = FALSE)
tmp <- tempfile()
file.create(paste0(tmp, c(".fst", ".log")))
regex <- "Weighted Fst estimate: (.*)"

################################################################################

test_that("snp_fst() works", {

  bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
  obj.bed <- bed(bedfile)

  pop <- rep(1:3, c(143, 167, 207))
  ind_pop <- split(seq_along(pop), pop)
  list_df_af <- lapply(ind_pop, function(ind) bed_MAF(obj.bed, ind.row = ind))

  lapply(list(1:3, 1:2, 2:3, c(1, 3)), function(pop_grps) {

    # PLINK to get Fst
    ind_pop_grp <- unlist(ind_pop[pop_grps])
    bigsnpr:::write.table2(cbind(obj.bed$fam[1:2], pop)[ind_pop_grp, ],
                           paste0(tmp, ".txt"))
    system(glue::glue("{plink} --bfile \"{sub_bed(bedfile)}\"",
                      " --fst --within {tmp}.txt --out {tmp}"),
           ignore.stdout = TRUE, ignore.stderr = TRUE)

    all_fst <- bigreadr::fread2(paste0(tmp, ".fst"))$FST
    fst <- as.double(
      sub(regex, "\\1", grep(regex, readLines(paste0(tmp, ".log")), value = TRUE)))

    # bigsnpr to get Fst
    expect_equal(snp_fst(list_df_af[pop_grps]),                 all_fst, tolerance = 1e-3)
    expect_equal(snp_fst(list_df_af[pop_grps], overall = TRUE), fst,     tolerance = 1e-5)
  })

})

################################################################################

test_that("snp_fst() works with missing values", {

  bedfile <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  obj.bed <- bed(bedfile)

  list_df_af <- lapply(list(1:100, 101:200),
                       function(ind) bed_MAF(obj.bed, ind.row = ind))

  # PLINK to get Fst
  bigsnpr:::write.table2(cbind(obj.bed$fam[1:2], rep(1:2, each = 100)),
                         paste0(tmp, ".txt"))
  system(glue::glue("{plink} --bfile \"{sub_bed(bedfile)}\"",
                    " --fst --within {tmp}.txt --out {tmp}"),
         ignore.stdout = TRUE, ignore.stderr = TRUE)

  all_fst <- bigreadr::fread2(paste0(tmp, ".fst"))$FST
  fst <- as.double(
    sub(regex, "\\1", grep(regex, readLines(paste0(tmp, ".log")), value = TRUE)))

  # bigsnpr to get Fst
  expect_equal(snp_fst(list_df_af),                 all_fst, tolerance = 2e-3)
  expect_equal(snp_fst(list_df_af, overall = TRUE), fst,     tolerance = 2e-4)

})

################################################################################
