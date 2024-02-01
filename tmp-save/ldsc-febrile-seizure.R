DIR <- "../datasets"

tgz <- runonce::download_file(
  "https://biobanks.dk/GWAS/Skotte2022_Febrile_Seizures.txt.gz", dir = DIR)

library(dplyr)
sumstats <- bigreadr::fread2(tgz) %>%
  rename(chr = "chromosome", a1 = "eff_allele", a0 = "alt_allele",
         beta = "BETA", beta_se = "SE") %>%
  as_tibble() %>%
  print()

info <- readRDS(runonce::download_file(
  "https://figshare.com/ndownloader/files/37802721",
  dir = DIR, fname = "map_hm3_plus.rds"))

library(bigsnpr)
info_snp <- snp_match(sumstats, info)
# 6,799,883 variants to be matched.
# 6 ambiguous SNPs have been removed.
# 992,627 variants have been matched; 0 were flipped and 472,955 were reversed.

(Neff <- quantile(8 / info_snp$beta_se^2, 0.999))  # 14992
sd_af <- with(info_snp, sqrt(2 * af_UKBB * (1 - af_UKBB)))
sd_ss <- 2 / sqrt(Neff * info_snp$beta_se^2)

library(ggplot2)
theme_set(theme_bw(14))
ggplot(data.frame(sd_af, sd_ss) %>% slice_sample(n = 50e3)) +
  geom_point(aes(sd_af, sd_ss)) +
  geom_abline(color = "red")

ldsc <- with(info_snp, {
  snp_ldsc(ld, ld_size = nrow(info), chi2 = (beta / beta_se)^2,
           sample_size = Neff, blocks = block_id, ncores = nb_cores())
})
#         int      int_se          h2       h2_se
# 0.997784156 0.008743117 0.495954617 0.049973969

ldsc[3:4] * coef_to_liab(K_pop = 0.036, K_gwas = 0.5)
#         h2      h2_se
# 0.38209955 0.03850157

