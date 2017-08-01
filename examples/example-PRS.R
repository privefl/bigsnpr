test <- snp_attachExtdata()
G <- test$genotypes
n <- nrow(G)
m <- ncol(G)

# get some phenotypes (LTM)
set.seed(5)
X <- attach.BM(G)
M <- 100 # number of causal variants
K <- 0.3 # prevalence
ind.causal <- sample(m, M)
betas <- rnorm(M, sd = sqrt(0.8 / M))
y.g <- drop(scale(X[, ind.causal]) %*% betas)
rm(X)
(v <- var(y.g))
y.e <- rnorm(n, sd = sqrt(1 - v))
(y01 <- as.numeric((y.g + y.e) > qnorm(1 - K)))

# PCA -> covariables
obj.svd <- snp_autoSVD(G, infos.chr = test$map$chromosome,
                       infos.pos = test$map$physical.position)

# train and test set
ind.train <- sort(sample(n, 400))
ind.test <- setdiff(rows_along(G), ind.train) # 117

# GWAS
gwas.train <- big_univLogReg(G, y01.train = y01[ind.train],
                             ind.train = ind.train,
                             covar.train = obj.svd$u[ind.train, ])
# clumping
ind.keep <- snp_clumping(G, infos.chr = test$map$chromosome,
                         ind.row = ind.train,
                         S = abs(gwas.train$score))
# -log10(p-values) and thresolding
summary(lpS <- -predict(gwas.train))
thrs <- seq(0, 3.5, by = 0.5)
nb.pred <- sapply(thrs, function(thr) sum(lpS[ind.keep] > thr))

# PRS
prs <- snp_PRS(G, betas = gwas.train$estim,
               ind.test = ind.test,
               ind.keep = ind.keep,
               lpS = lpS,
               thr.list = thrs)

# AUC as a function of the number of predictors
aucs <- apply(prs, 2, AUC, target = y01[ind.test])
library(ggplot2)
bigstatsr:::MY_THEME(qplot(nb.pred, aucs)) +
  geom_line() +
  scale_x_log10(breaks = nb.pred) +
  labs(x = "Number of predictors", y = "AUC")
