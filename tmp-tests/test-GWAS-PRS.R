require(bigsnpr)

test <- snp_attach("backingfiles/celiac300_sub1.rds")
G <- test$genotypes
G@code <- c(0, 1, rep(2, 254))

ind.train <- sort(sample(rows_along(G), size = 10e3))
ind.test <- setdiff(rows_along(G), ind.train)

ind.excl <- snp_indLRLDR(test$map$chromosome, test$map$physical.pos)
ind.keep <- snp_clumping(G, test$map$chromosome, thr.r2 = 0.2,
                         ind.row = ind.train, exclude = ind.excl)
svd <- big_randomSVD(G, snp_scaleBinom(), ind.col = ind.keep, ncores = 3)

plot(svd$u)
plot(svd$u[, 3:4])
plot(svd$d)

y01 <- bigstatsr:::transform_levels(test$fam$affection)
gwas <- big_univLogReg(G, y01.train = y01[ind.train], ind.train = ind.train,
                       covar.train = svd$u[ind.train, 1:2], ncores2 = 2)

snp_manhattan(gwas, test$map)

snp_manhattan(bigsnpr:::snp_gc(gwas), test$map)
snp_qq(gwas)
snp_qq(snp_gc(gwas))

ind.keep2 <- snp_clumping(G, test$map$chromosome, ind.train,
                          S = abs(gwas$score))
prs <- snp_PRS(G, gwas$estim, ind.test, ind.keep = ind.keep2,
               lpS = -log10(predict(gwas)), thr.list = 0:20 / 2)
aucs <- apply(prs, 2, AUC, target = y01[ind.test])
plot(aucs)

prs2 <- snp_PRS(G, gwas$estim, ind.test, ind.keep = ind.keep2,
               lpS = -log10(predict(snp_gc(gwas))), thr.list = 0:20 / 2)
aucs2 <- apply(prs2, 2, AUC, target = y01[ind.test])
plot(aucs2)


ind.keep3 <- snp_clumping(G, test$map$chromosome, ind.train,
                          S = abs(gwas$score), thr.r2 = 0.2)
prs3 <- snp_PRS(G, gwas$estim, ind.test, ind.keep = ind.keep3,
                lpS = -log10(predict(gwas)), thr.list = 0:20 / 2)
aucs3 <- apply(prs3, 2, AUC, target = y01[ind.test])
plot(aucs3)

test2 <- big_CMSA(big_spSVM, feval = AUC, y.train = y01[ind.train],
                 ind.train = ind.train, covar.train = svd$u[ind.train, 1:2],
                 X. = G, K = 10, ncores = 3, alpha = 0.5)
svm <- big_prodVec(G, test2[cols_along(G)], ind.row = ind.test)
aucs4 <- AUC(svm, target = y01[ind.test])
svm2 <- svd$u[ind.test, 1:2] %*% tail(test, 2)
AUC(svm2, target = y01[ind.test])
AUC(svm + svm2, target = y01[ind.test])

require(ggplot2)
pop.files <- paste0("../Dubois2010_data/others/cluster_",
                    c("IT", "NL", "UK1_v2", "UK3_v2", "Finn_v2"))
test <- snp_getPops(test, pop.files = pop.files,
                    col.sample.ID = 2, col.family.ID = 3)
POPS <- c("Netherlands", "Italy", "UK1", "UK2", "Finland")
ggplot(data = data.frame(score = svm + svm2,
                         status = ifelse(y01[ind.test], "case", "control"))) +
  geom_density(aes(x = score, group = status, fill = status), alpha = 0.3) +
  facet_wrap(~ POPS[test$fam$family.ID[ind.test]])
