sumstats <- bigreadr::fread2("tests/test_data/sim1_0_ss.txt")
hist(sumstats$PVAL)

test.bed <- "tests/test_data/sim1_0_test.bed"
library(bigsnpr)
snp_readBed(test.bed)
obj.bigsnp <- snp_attach("tests/test_data/sim1_0_test.rds")
G <- obj.bigsnp$genotypes
dim(G) # 2000 2000

K <- big_cor(G)
diag0 <- diag(K)

m <- ncol(G)
n <- nrow(G)
h2 <- 0.1
diag(K) <- diag0 + m / (n * h2)
new_beta <- solve(K[], sumstats$BETA)
plot(sumstats$BETA, new_beta, pch = 20); abline(0, 1, col = "red")

pred <- big_prodVec(G, new_beta)
y <- obj.bigsnp$fam$affection
cor(-pred, y)^2

CHR <- obj.bigsnp$map$chromosome
POS <- obj.bigsnp$map$physical.pos
lpval <- -log10(sumstats$PVAL)
NCORES <- nb_cores()
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
