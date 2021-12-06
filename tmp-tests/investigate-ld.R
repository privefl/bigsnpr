
library(dplyr)
library(bigsnpr)

sumstats <- bigreadr::fread2(
  "tmp-data/Mahajan.NatGenet2018b.T2D.European/Mahajan.NatGenet2018b.T2D.European.txt",
  select = c("Chr", "Pos", "EA", "NEA", "EAF", "Beta", "SE", "Neff"),
  col.names = c("chr", "pos", "a1", "a0", "freq", "beta", "beta_se", "n_eff")
) %>%
  filter(n_eff > 230e3, pmin(freq, 1 - freq) > 0.01)

info <- readRDS(runonce::download_file(
  "https://ndownloader.figshare.com/files/25503788",
  dir = "tmp-data", fname = "map_hm3_ldpred2.rds")) %>%
  filter(chr == 8)

info_snp <- snp_match(sumstats, info)
ind <- info_snp$`_NUM_ID_`

corr0 <- readRDS(
  runonce::download_file("https://figshare.com/ndownloader/files/24927941",
                         dir = "tmp-data", fname = "LD_chr8.rds"))

corr <- as_SFBM(corr0[ind, ind])

auto <- snp_ldpred2_auto(corr, info_snp, h2_init = 0.05, vec_p_init = 0.01,
                         verbose = TRUE, burn_in = 100, num_iter = 100)[[1]]
# seems ok
