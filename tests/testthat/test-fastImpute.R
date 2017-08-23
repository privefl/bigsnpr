################################################################################

context("FAST_IMPUTE")

################################################################################

bigsnp <- snp_attachExtdata()
G <- bigsnp$genotypes
expect_equal(G@code, bigsnpr:::CODE_012)
corr <- snp_cor(G)
library(Matrix)
ind <- which(apply(corr, 2, function(x) max(x[x < 1])) > 0.6)

library(foreach)
indNA <- foreach(j = ind, .combine = "rbind") %do% {
  cbind(sample(nrow(G), 100), j)
}

# Create a copy of our bigSNP example
library(magrittr)
tmpfile <- tempfile()
bigsnp.copy <- snp_writeBed(bigsnp, bedfile = paste0(tmpfile, ".bed")) %>%
  snp_readBed(backingfile = basename(tmpfile),
              backingpath = dirname(tmpfile)) %>%
  snp_attach()
# Fill some missing values
GNA <- bigsnp.copy$genotypes
XNA <- attach.BM(GNA)
XNA[indNA] <- as.raw(3)

################################################################################

# Fast imputation
time <- system.time(
  infos <- snp_fastImpute(GNA, infos.chr = bigsnp$map$chromosome,
                          p.train = 0.6)
)
# expect_lt(time[3], 10)

XNA@code <- bigsnpr:::CODE_IMPUTE_PRED
library(dplyr)
infosNA <- data_frame(
  col = indNA[, 2],
  error =  XNA[indNA] != attach.BM(G)[indNA]
) %>%
  group_by(col) %>%
  summarise(nb_err = sum(error)) %>%
  mutate(nb_err_est = (infos$pError * infos$pNA * nrow(G))[col])

pval <- anova(lm(nb_err ~ nb_err_est - 1, data = infosNA))$`Pr(>F)`[[1]]
expect_lt(pval, 2e-16)

################################################################################
