library(bigsnpr)
library(bigreadr)
library(dplyr)

NCORES <- 4  # TO MODIFY

## Information for the variants provided in the LD reference
map_ldref <- readRDS(url("https://figshare.com/ndownloader/files/25503788"))

# https://drive.google.com/file/d/1lmnpNDpPiwliMvd2Bb33hPpzzVJkLFp9/view?usp=share_link
sumstats <- fread2("../datasets/GCST90086597_buildGRCh37.tsv.gz",
                   select = c("chromosome", "base_pair_location",
                              "other_allele", "effect_allele",
                              "beta", "standard_error", "effect_allele_frequency"),
                   col.names = c("chr", "pos", "a0", "a1", "beta", "beta_se", "freq")) %>%
  filter(chr %in% 1:22) %>%
  mutate(chr = as.integer(chr))

sumstats$n_eff <- 5368


# Ancestry inference
all_freq <- bigreadr::fread2(
  runonce::download_file("https://figshare.com/ndownloader/files/31620968",
                         dir = "../datasets", fname = "ref_freqs.csv.gz"))
projection <- bigreadr::fread2(
  runonce::download_file("https://figshare.com/ndownloader/files/31620953",
                         dir = "../datasets", fname = "projection.csv.gz"))
matched <- snp_match(sumstats, all_freq[1:5], return_flip_and_rev = TRUE) %>%
  mutate(freq = ifelse(`_REV_`, 1 - freq, freq))
# 7,311,951 variants to be matched.
# 32 ambiguous SNPs have been removed.
# 4,349,781 variants have been matched; 0 were flipped and 1,116,879 were reversed.

res <- snp_ancestry_summary(
  freq = matched$freq,
  info_freq_ref = all_freq[matched$`_NUM_ID_`, -(1:5)],
  projection = projection[matched$`_NUM_ID_`, -(1:5)],
  correction = c(1, 1, 1, 1.008, 1.021, 1.034, 1.052, 1.074, 1.099,
                 1.123, 1.15, 1.195, 1.256, 1.321, 1.382, 1.443)
)
# Note that some ancestry groups from the reference are very close to one another,
# and should be merged a posteriori.
group <- colnames(all_freq)[-(1:5)]
group[group %in% c("Scandinavia", "United Kingdom", "Ireland")]   <- "Europe (North West)"
group[group %in% c("Europe (South East)", "Europe (North East)")] <- "Europe (East)"
grp_fct <- factor(group, unique(group))
final_res <- tapply(res, grp_fct, sum)
round(100 * final_res[final_res > 0.001], 1)
# Finland Europe (North West)
#     1.7                98.3


info_snp <- snp_match(sumstats, map_ldref)
# 7,311,951 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 1,034,011 variants have been matched; 0 were flipped and 313,434 were reversed.
(info_snp <- tidyr::drop_na(tibble::as_tibble(info_snp)))

(sd_y <- with(info_snp, sqrt(0.5 * quantile(n_eff * beta_se^2 + beta^2, 0.01))))
# 0.9843047

sd_ldref <- with(info_snp, sqrt(2 * af_UKBB * (1 - af_UKBB)))
sd_af <- with(info_snp, sqrt(2 * freq * (1 - freq)))
sd_ss <- with(info_snp, 1 / sqrt(n_eff * beta_se^2 + beta^2))

is_bad <-
  sd_ss < (0.7 * sd_af) | sd_ss > (sd_af + 0.1) | sd_ss < 0.1 | sd_af < 0.05

library(ggplot2)
qplot(sd_af, sd_ss, color = is_bad, alpha = I(0.5),
      data = slice_sample(data.frame(sd_af, sd_ss, is_bad), n = 50e3)) +
  theme_bigstatsr() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations derived from allele frequencies",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")
qplot(sd_ldref, sd_ss, color = is_bad, alpha = I(0.5),
      data = slice_sample(data.frame(sd_ldref, sd_ss, is_bad), n = 50e3)) +
  theme_bigstatsr() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations derived from allele frequencies of the LD reference",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")

df_beta <- info_snp[!is_bad, ]


tmp <- tempfile(tmpdir = "tmp-data")

for (chr in 1:22) {

  cat(chr, ".. ", sep = "")

  ## indices in 'df_beta'
  ind.chr <- which(df_beta$chr == chr)
  ## indices in 'map_ldref'
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  ## indices in 'corr_chr'
  ind.chr3 <- match(ind.chr2, which(map_ldref$chr == chr))

  corr_chr <- readRDS(paste0("../datasets/ldref/LD_with_blocks_chr", chr, ".rds"))[ind.chr3, ind.chr3]

  if (chr == 1) {
    corr <- as_SFBM(corr_chr, tmp, compact = TRUE)
  } else {
    corr$add_columns(corr_chr, nrow(corr))
  }
}

# Heritability estimation of LD score regression
# to be used as a starting value in LDpred2-auto
(ldsc <- with(df_beta, snp_ldsc(ld, ld_size = nrow(map_ldref),
                                chi2 = (beta / beta_se)^2,
                                sample_size = n_eff,
                                ncores = NCORES)))
#        int     int_se         h2      h2_se
# 1.17693343 0.01450483 2.83327279 2.82301985


# LDpred2-auto
coef_shrink <- 0.95
multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = ldsc[["h2"]],
                               vec_p_init = seq_log(1e-4, 0.2, length.out = 20),
                               burn_in = 500, num_iter = 500, report_step = 20,
                               allow_jump_sign = FALSE, shrink_corr = coef_shrink,
                               ncores = NCORES)
(range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))

#  [1] 1.604992 1.609106 1.607678 1.606848 1.604719 1.605809 1.602778 1.608586
#  [9] 1.602154 1.596995 1.599681 1.612114 1.607826 1.602184 1.608860 1.603382
# [17] 1.600513 1.607723 1.605344 1.61388

(keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE)) & range < 3))

all_h2 <- sapply(multi_auto[keep], function(auto) tail(auto$path_h2_est, 500))
quantile(all_h2, c(0.5, 0.025, 0.975))

#      50%     2.5%    97.5%
# 3.090475 2.578516 3.580651

all_p <- sapply(multi_auto[keep], function(auto) tail(auto$path_p_est, 500))
quantile(all_p, c(0.5, 0.025, 0.975))

#         50%        2.5%       97.5%
# 0.005265126 0.002900875 0.008655001

all_alpha <- sapply(multi_auto[keep], function(auto) tail(auto$path_alpha_est, 500))
quantile(all_alpha, c(0.5, 0.025, 0.975))
#        50%       2.5%     97.5%
# -1.0012462 -1.1300291 -0.862035

# cleanup
rm(corr); gc(); file.remove(paste0(tmp, ".sbk"))

r2 <- with(sumstats, beta^2 / (beta^2 + n_eff * beta_se^2))
ind <- which(r2 > 0.01)
plot(ind, r2[ind])
sumstats[ind, ]

r2 <- with(df_beta, beta^2 / (beta^2 + n_eff * beta_se^2))
ind <- which(r2 > 0.6)
plot(ind, r2[ind])
df_beta[ind, ]
df_beta$beta[ind]

round(100 * corr[ind, ind], 1)
