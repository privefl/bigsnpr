test <- snp_attachExtdata()
G <- test$genotypes
y01 <- test$fam$affection - 1

# PCA -> covariables
obj.svd <- snp_autoSVD(G, infos.chr = test$map$chromosome,
                       infos.pos = test$map$physical.position)

# train and test set
ind.train <- sort(sample(nrow(G), 400))
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
summary(lpS.keep <- -predict(gwas.train)[ind.keep])
thrs <- seq(0, 3.5, by = 0.5)
nb.pred <- sapply(thrs, function(thr) sum(lpS.keep > thr))

# PRS
prs <- snp_PRS(G, betas.keep = gwas.train$estim[ind.keep],
               ind.test = ind.test,
               ind.keep = ind.keep,
               lpS.keep = lpS.keep,
               thr.list = thrs)

# AUC as a function of the number of predictors
aucs <- apply(prs, 2, AUC, target = y01[ind.test])
library(ggplot2)
bigstatsr:::MY_THEME(qplot(nb.pred, aucs)) +
  geom_line() +
  scale_x_log10(breaks = nb.pred) +
  labs(x = "Number of predictors", y = "AUC")
