require(bigsnpr)

test <- snp_attach("backingfiles/popresNA.rds")
G <- test$genotypes
G@code <- bigsnpr:::CODE_DOSAGE

ind.keep <- snp_pruning(G, test$map$chromosome, thr.r2 = 0.1)
G@code <- snp_pruning()


pca <- big_randomSVD(G, snp_scaleBinom(), ind.col = ind.keep, k = 20)

require(SNPRelate)

genofile <- snpgdsOpen("../../Bureau/POPRES_data/popresSub.gds")

print(system.time(
  relate <- snpgdsIBDMoM(genofile, num.thread = 1,
                         snp.id = test$map$marker.ID[ind.keep])
))
pi_hat <- (1 - relate$k0 - relate$k1 / 2)



ms <- snp_scaleBinom()(G, ind.col = ind.keep)
X <- attach.BM(G)[, ind.keep] # 1 Gb
X <- bigstatsr:::scaling(X, ms$mean, ms$sd)
mean(X[, 1])
sd(X[, 1])
print(system.time(
  xxt <- tcrossprod(X)
))
Rcpp::sourceCpp('tmp-save/center.cpp')
diags <- sqrt(diag(xxt))
xxt <- toCorr(xxt, sqrt_diags = diags)
xxt[1:5, 1:5]
all.kinship <- xxt[upper.tri(xxt)]

round(xxt[1:8, 1:8], 3)
round(pi_hat[1:8, 1:8], 3)
ind <- which(xxt > 0.06, arr.ind = TRUE)
ind2 <- ind[ind[, 1] < ind[, 2], ]
plot(pi_hat[ind2], xxt[ind2]); abline(0, 1, col = "red")



### small number of indices
keep.sample <- sort(sample(length(ind.keep), 10e3))
X2 <- attach.BM(G)[, ind.keep[keep.sample]]
X2 <- bigstatsr:::scaling(X2, ms$mean[keep.sample], ms$sd[keep.sample])
mean(X2[, 1])
sd(X2[, 1])
print(system.time(
  xxt2 <- tcrossprod(X2)
))
diags2 <- sqrt(diag(xxt2))
xxt2 <- toCorr(xxt2, sqrt_diags = diags2)
ind <- which(xxt2 > 0.04, arr.ind = TRUE)
ind2 <- ind[ind[, 1] < ind[, 2], ]
plot(xxt2[ind2], xxt[ind2])


# first 20
U <- pca$u
plot(U)
tmp <- dist(U)
round(as.matrix(tmp)[1:10, 1:10], 3)
round(pi_hat[1:10, 1:10], 3)
ind <- which(pi_hat > 0.05, arr.ind = TRUE)
ind2 <- ind[ind[, 1] < ind[, 2], ]
plot(pi_hat[ind2], as.matrix(tmp)[ind2])

# first 3
U <- pca$u[, c(1:5)]
plot(U)
tmp <- dist(U)
round(as.matrix(tmp)[1:10, 1:10], 3)
round(pi_hat[1:10, 1:10], 3)
ind <- which(pi_hat > 0.05, arr.ind = TRUE)
ind2 <- ind[ind[, 1] < ind[, 2], ]
plot(pi_hat[ind2], as.matrix(tmp)[ind2])
mean(as.matrix(tmp) > 0.04) # 0.76
