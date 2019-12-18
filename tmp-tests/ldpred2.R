sumstats <- bigreadr::fread2("tests/test_data/sim1_0_ss.txt")
hist(sumstats$PVAL)

test.bed <- "tests/test_data/sim1_0_test.bed"
library(bigsnpr)
snp_readBed(test.bed)
obj.bigsnp <- snp_attach("tests/test_data/sim1_0_test.rds")
G <- obj.bigsnp$genotypes
dim(G) # 2000 2000
CHR <- obj.bigsnp$map$chromosome
POS <- 1000 * obj.bigsnp$map$physical.pos + 1
lpval <- -log10(sumstats$PVAL)
NCORES <- nb_cores()

K <- snp_cor(G, infos.pos = POS, alpha = 1)

ld <- apply(K, 2, crossprod)
all.equal(qchisq(sumstats$PVAL, df = 1, lower.tail = FALSE),
          qchisq(-lpval * log(10), df = 1, lower.tail = FALSE, log.p = TRUE))
chi2 <- qchisq(-lpval * log(10), df = 1, lower.tail = FALSE, log.p = TRUE)

library(dplyr)
data.frame(ld, chi2) %>%
  group_by(cut(ld, c(-Inf, quantile(ld, 1:19 / 20), Inf))) %>%
  summarise_all(mean) %>%
  plot(chi2 ~ ld, data = .)
(ldsc <- LDSC(chi2 = chi2, ld = ld, M = ncol(G), N = sumstats$N[[1]]))

diag(K) <- 1 + 1 / ldsc[["SLP"]]
new_beta <- as.vector(Matrix::solve(K, sumstats$BETA))
plot(sumstats$BETA, new_beta, pch = 20); abline(0, 1, col = "red")

pred <- big_prodVec(G, new_beta)
y <- obj.bigsnp$fam$affection
cor(-pred, y)^2

corr <- sapply(10^(-6:0), function(alpha) {
  K <- snp_cor(G, infos.pos = POS, alpha = alpha, size = 10e3)

  ld <- apply(K, 2, crossprod)
  all.equal(qchisq(sumstats$PVAL, df = 1, lower.tail = FALSE),
            qchisq(-lpval * log(10), df = 1, lower.tail = FALSE, log.p = TRUE))
  chi2 <- qchisq(-lpval * log(10), df = 1, lower.tail = FALSE, log.p = TRUE)

  (ldsc <- LDSC(chi2 = chi2, ld = ld, M = ncol(G), N = sumstats$N[[1]]))

  diag(K) <- 1 + 1 / ldsc[["SLP"]]
  new_beta <- as.vector(Matrix::solve(K, sumstats$BETA))
  plot(sumstats$BETA, new_beta, pch = 20, main = alpha); abline(0, 1, col = "red")

  pred <- big_prodVec(G, new_beta)
  y <- obj.bigsnp$fam$affection
  cor(-pred, y)^2
})
corr
# size=10e3: 0.06875931 0.06912793 0.06894654 0.06933440 0.06895514 0.06807963 0.06673039
# size=500:  0.06875931 0.06910886 0.06886562 0.06905332 0.06852254 0.06819528 0.06918475
# size=100:  0.06875931 0.06911954 0.06887603 0.06916707 0.06900139 0.06886611 0.06874929
# size=10:   0.06538603 0.06538603 0.06538603 0.06538603 0.06537984 0.06535629 0.06530207

# Clumping
all_keep <- snp_grid_clumping(G, CHR, 1000 * POS + 1, lpS = lpval, ncores = NCORES,
                              grid.base.size = 200)
attr(all_keep, "grid")
# Thresholding
beta <- -sumstats$BETA
tmp <- tempfile()
multi_PRS <- snp_grid_PRS(G, all_keep, beta, lpval,
                          backingfile = tmp,
                          n_thr_lpS = 50, ncores = NCORES)
dim(multi_PRS)  ## 2000 x 350

library(tidyverse)
grid2 <- attr(all_keep, "grid") %>%
  mutate(thr.lp = list(attr(multi_PRS, "grid.lpS.thr")), num = row_number()) %>%
  unnest()

## Warning: `cols` is now required.
## Please use `cols = c(thr.lp)`

s <- nrow(grid2)
grid2$r2 <- big_apply(multi_PRS, a.FUN = function(X, ind, s, y.train) {
  # Sum over all chromosomes, for the same C+T parameters
  print(single_PRS <- rowSums(X[, ind, drop = FALSE]))  ## replace by 0:21 in real data
  cor(single_PRS, y.train)^2
}, ind = 1:s, s = s, y.train = y,
a.combine = 'c', block.size = 1, ncores = 1)

max_prs <- grid2 %>% arrange(desc(r2)) %>% slice(1:10) %>% print() %>% slice(1)

K2 <- snp_cor(G, infos.pos = 1000 * POS + 1)
str(K2)
K2[1:5, 1:8]
K[1:5, 1:8]
K2[, 1:2]
