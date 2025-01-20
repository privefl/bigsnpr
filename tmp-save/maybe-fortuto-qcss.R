library(testthat)
library(bigsnpr)

bedfile <- file.path(tempdir(), "tmp-data/public-data3.bed")
if (!file.exists(rdsfile <- sub_bed(bedfile, ".rds"))) {
  zip <- tempfile(fileext = ".zip")
  download.file(
    "https://github.com/privefl/bigsnpr/blob/master/data-raw/public-data3.zip?raw=true",
    destfile = zip, mode = "wb")
  unzip(zip, exdir = tempdir())
  rds <- snp_readBed(bedfile)
  expect_identical(normalizePath(rds), normalizePath(rdsfile))
}

obj.bigSNP <- snp_attach(rdsfile)
G <- obj.bigSNP$genotypes
y <- obj.bigSNP$fam$affection
POS2 <- obj.bigSNP$map$genetic.dist + 1000 * obj.bigSNP$map$chromosome

sumstats <- bigreadr::fread2(file.path(tempdir(), "tmp-data/public-data3-sumstats.txt"))
sumstats$n_eff <- sumstats$N
map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
df_beta <- snp_match(sumstats, map, join_by_pos = FALSE)

ind_var <- df_beta$`_NUM_ID_`
corr0 <- snp_cor_extendedThr(
  G, ind.col = ind_var, thr_r2 = 0.2,
  infos.pos = POS2[ind_var], size = 1 / 1000, ncores = 2)
corr <- as_SFBM(corr0)

maf <- snp_MAF(G, ind.col = ind_var)
sd0 <- sqrt(2 * maf * (1 - maf))

expect_equal(
  snp_qcimp_sumstats(corr, df_beta, sd0 = sd0, max_run = 2),
  snp_qcimp_sumstats(corr, df_beta, sd0 = sd0, max_run = 2)
)
#
# expect_equal(
#   snp_qcimp_sumstats(corr, df_beta, sd0 = sd0, max_run = 1),
#   snp_qcimp_sumstats(corr, df_beta, sd0 = sd0, max_run = 1, ncores = 2)
# )
# 1 iter is fine
# 2 iter is fine with ncores = 1
# any is fine with ncores = 1
# the problem comes from ncores = 2 and max_run > 1

# test0 <- replicate(10, snp_qcimp_sumstats(corr, df_beta, sd0 = sd0, max_run = 2))
# #18378 / 10129 / 19061 is TRUE

expect_equal(
  test1 <- snp_qcimp_sumstats(corr, df_beta, sd0 = sd0, max_run = 2, ncores = 3),
  test2 <- snp_qcimp_sumstats(corr, df_beta, sd0 = sd0, max_run = 2, ncores = 3)
)

# only r2_max is the same
# thr are all 0.2 / keep are all TRUE

expect_equal(
  test3 <- snp_qcimp_sumstats(corr, df_beta, sd0 = sd0, ncores = 3),
  test4 <- snp_qcimp_sumstats(corr, df_beta, sd0 = sd0, ncores = 3)
)

plot(lengths(test1$ind_imp), lengths(test2$ind_imp))

res <- snp_qcimp_sumstats(corr, df_beta, sd0 = sd0, ncores = 2)
expect_false(any(res$rm_qc, na.rm = TRUE))
expect_equal(snp_qcimp_sumstats(corr, df_beta,  sd0 = sd0, ncores = 2), res)
