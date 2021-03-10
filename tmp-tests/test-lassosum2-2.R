# Load packages bigsnpr and bigstatsr
library(bigsnpr)
# Read from bed/bim/fam, it generates .bk and .rds files.
# snp_readBed("tmp-data/public-data.bed")
# Attach the "bigSNP" object in R session
obj.bigSNP <- snp_attach("tmp-data/public-data.rds")
# See how the file looks like
str(obj.bigSNP, max.level = 2, strict.width = "cut")
# Get aliases for useful slots
G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
y   <- obj.bigSNP$fam$affection - 1
NCORES <- nb_cores()
# Read external summary statistics
sumstats <- bigreadr::fread2("tmp-data/public-data-sumstats.txt")
str(sumstats)


set.seed(1)
ind.val <- sample(nrow(G), 400)
ind.test <- setdiff(rows_along(G), ind.val)

sumstats$n_eff <- 4 / (1 / sumstats$n_case + 1 / sumstats$n_control)
sumstats$n_case <- sumstats$n_control <- NULL
names(sumstats) <- c("chr", "rsid", "pos", "a0", "a1", "beta", "beta_se", "p", "n_eff")
map <- obj.bigSNP$map[-(2:3)]
names(map) <- c("chr", "pos", "a0", "a1")
info_snp <- snp_match(sumstats, map, strand_flip = FALSE)

POS2 <- snp_asGeneticPos(CHR, POS, dir = "tmp-data", ncores = NCORES)

ind.chr <- which(info_snp$chr == 2)
df_beta <- info_snp[ind.chr, c("beta", "beta_se", "n_eff")]
ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
corr0 <- snp_cor(G, ind.col = ind.chr2, ncores = NCORES,
                 infos.pos = POS2[ind.chr2], size = 3 / 1000)
tmp <- tempfile(tmpdir = "tmp-data")
corr <- bigsparser::as_SFBM(as(corr0, "dgCMatrix"), tmp)

# takes several minutes if you do not have many cores
beta_grid <- snp_lassosum2(corr, df_beta, ncores = NCORES)

pred_grid <- big_prodMat(G, beta_grid, ind.col = ind.chr2)

lambda0 <- max(abs(df_beta$beta / sqrt(df_beta$n_eff * df_beta$beta_se^2 + df_beta$beta^2)))
params <- expand.grid(
  lambda = seq_log(0.01 * lambda0, lambda0, 20),
  s = c(0.2, 0.5, 0.8, 0.9, 0.95, 1)
)
params$score <- big_univLogReg(as_FBM(pred_grid[ind.val, ]), y[ind.val])$score

library(ggplot2)
ggplot(params, aes(x = lambda, y = score, color = as.factor(s))) +
  theme_bigstatsr() +
  geom_point() +
  geom_line() +
  scale_x_log10(breaks = 10^(-5:0)) +
  labs(y = "GLM Z-Score", color = "s") +
  theme(legend.position = "top", panel.spacing = unit(1, "lines"))


params$sparsity <- colMeans(beta_grid == 0)

ggplot(params, aes(x = lambda, y = sparsity, color = as.factor(s))) +
  theme_bigstatsr() +
  geom_point() +
  geom_line() +
  scale_x_log10(breaks = 10^(-5:0)) +
  labs(y = "Sparsity", color = "s") +
  theme(legend.position = "top", panel.spacing = unit(1, "lines"))
