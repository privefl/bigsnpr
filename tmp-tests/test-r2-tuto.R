## ------------------------------------------------------------------------
# Load packages bigsnpr and bigstatsr
library(bigsnpr)
# Read from bed/bim/fam, it generates .bk and .rds files.
if (!file.exists("tmp-data/public-data3.bk"))
  snp_readBed("tmp-data/public-data3.bed")
# Attach the "bigSNP" object in R session
obj.bigSNP <- snp_attach("tmp-data/public-data3.rds")
# See how the file looks like
str(obj.bigSNP, max.level = 2, strict.width = "cut")
# Get aliases for useful slots
G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
y   <- obj.bigSNP$fam$affection
(NCORES <- nb_cores())


## ------------------------------------------------------------------------
# Read external summary statistics
sumstats <- bigreadr::fread2("tmp-data/public-data3-sumstats.txt")
str(sumstats)


## ------------------------------------------------------------------------
set.seed(1)
ind.val <- sample(nrow(G), 350)
ind.test <- setdiff(rows_along(G), ind.val)


## ---- error=TRUE---------------------------------------------------------
# sumstats$n_eff <- 4 / (1 / sumstats$n_case + 1 / sumstats$n_control)
# sumstats$n_case <- sumstats$n_control <- NULL
sumstats$n_eff <- sumstats$N
map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
df_beta <- snp_match(sumstats, map, join_by_pos = FALSE)  # use rsid instead of pos


## ------------------------------------------------------------------------
# To convert physical positions (in bp) to genetic positions (in cM), use
# POS2 <- snp_asGeneticPos(CHR, POS, dir = "tmp-data", ncores = NCORES)
# To avoid downloading "large" files, `POS2` has been precomputed here
POS2 <- obj.bigSNP$map$genetic.dist


## ------------------------------------------------------------------------
tmp <- tempfile(tmpdir = "tmp-data")

for (chr in 1:22) {

  # print(chr)

  ## indices in 'df_beta'
  ind.chr <- which(df_beta$chr == chr)
  ## indices in 'G'
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]

  corr0 <- snp_cor(G, ind.col = ind.chr2, size = 3 / 1000,
                   infos.pos = POS2[ind.chr2], ncores = NCORES)

  if (chr == 1) {
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, tmp, compact = TRUE)
  } else {
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
}


## ------------------------------------------------------------------------
file.size(corr$sbk) / 1024^3  # file size in GB


## ------------------------------------------------------------------------
# Estimate of h2 from LD Score regression
(ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                sample_size = n_eff, blocks = NULL)))
h2_est <- ldsc[["h2"]]



## ---- cache=TRUE---------------------------------------------------------
coef_shrink <- 0.95
multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                               vec_p_init = seq_log(1e-4, 0.2, length.out = 50),
                               burn_in = 500, num_iter = 500, report_step = 20,
                               allow_jump_sign = FALSE, shrink_corr = coef_shrink,
                               ncores = NCORES)
(range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
(keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE)) & range < 3))


## ------------------------------------------------------------------------
bsamp <- lapply(multi_auto[keep], function(auto) auto$sample_beta)
all_r2 <- do.call("cbind", lapply(seq_along(bsamp), function(ic) {
  b1 <- bsamp[[ic]]
  Rb1 <- apply(b1, 2, function(x)
    coef_shrink * bigsparser::sp_prodVec(corr, x) + (1 - coef_shrink) * x)
  b2 <- do.call("cbind", bsamp[-ic])
  b2Rb1 <- as.matrix(Matrix::crossprod(b2, Rb1))
}))
quantile(all_r2, c(0.5, 0.025, 0.975))


## ------------------------------------------------------------------------
beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
pred_auto <- big_prodVec(G, beta_auto, ind.col = df_beta[["_NUM_ID_"]])
pcor(pred_auto, y, NULL)^2


scale <- with(df_beta, sqrt(n_eff * beta_se^2 + beta^2))
beta <- beta_auto / scale
bRb <- crossprod(beta, bigsparser::sp_prodVec(corr, beta))[1]
quantile(all_r2^2 / bRb, c(0.5, 0.025, 0.975))


## ------------------------------------------------------------------------
postp <- rowMeans(sapply(multi_auto[keep], function(auto) auto$postp_est))
library(ggplot2); qplot(y = postp, alpha = I(0.2)) + theme_bigstatsr()





## ------------------------------------------------------------------------
# Some cleaning
rm(corr); gc(); file.remove(paste0(tmp, ".sbk"))

