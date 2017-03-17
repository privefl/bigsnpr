require(bigsnpr)

popres <- snp_attach("backingfiles/popresNA_sub1.rds")
X <- attach.BM(popres$genotypes)
n <- nrow(X)
m <- ncol(X)

# print(table(
#   nbNA <- VGAM::rbetabinom.ab(m, size = n, shape1 = 0.6, shape2 = 100)
# ))
#
# for (j in 1:m) {
#   indNA <- sample(n, size = nbNA[j])
#   X[indNA, j] <- as.raw(X[indNA, j] + 4)
# }
X[, 1]
#
X2 <- X
X2@code <- bigsnpr:::CODE_DOSAGE
X2[, 1]
print(store.NA <- readRDS("storeNA_popres.rds")) #X2[which(is.na(X[,]))]

# correlation between SNPs
Rcpp::sourceCpp('src/corr.cpp')
q.alpha <- stats::qchisq(0.2, df = 1, lower.tail = FALSE)
require(Matrix)
m.part <- 1000
corr <- symmpart(corMat(X, rowInd = 1:n, colInd = 1:1000,
                        size = 1000, thr = q.alpha / 1:n))

#### IMPUTATION TEST ####

require(xgboost)

nbNA <- integer(m.part)
error <- rep(NA_real_, m.part)
error.save <- rep(2, m.part)
indNA.part <- which(is.na(X[, 1:m.part]))

# useful functions
round2 <- function(pred) as.raw(round(100 * pmin(pmax(pred, 0), 2)) + 7)
round3 <- function(pred) (pred > 0.5) + (pred > 1.5)
round4 <- function(pred) round2(0:2 %*% matrix(pred, 3))
round5 <- function(pred) apply(matrix(pred, 3), 2, which.max) - 1

# X2.svd <- big_SVD(X2, snp_scaleBinom(), k = 6)

# arr.indNA <- which(is.na(X[, 1:1000]), arr.ind = TRUE)
# for (i in 1:nrow(arr.indNA))
#   X[arr.indNA[i, 1], arr.indNA[i, 2]] <- as.raw(3)
# X2[, 1]


# imputation
for (i in 1:m.part) {
  if (!(i %% 10)) print(i)
  X.label <- X[, i]
  nbNA[i] <- l <- length(indNA <- which(is.na(X.label)))
  if (l > 0) {
    indNoNA <- setdiff(1:n, indNA)
    ind.train <- sort(sample(indNoNA, 0.8 * length(indNoNA)))
    ind.val <- setdiff(indNoNA, ind.train)

    X.data <- X2[, which(corr[, i] != 0), drop = FALSE]

    bst <- xgboost(data = X.data[ind.train, , drop = FALSE],
                   label = X.label[ind.train],
                   objective = "multi:softprob",
                   num_class = 3,
                   base_score = mean(X.label[ind.train]),
                   nrounds = round(runif(1, 40, 50)),
                   params = list(max_depth = sample(3:4, 1), gamma = 1,
                                 colsample_bytree = runif(1, 0.7, 1),
                                 colsample_bylevel = runif(1, 0.6, 1)),
                   nthread = 1,
                   verbose = 0,
                   save_period = NULL)

    # error of validation
    pred2 <- predict(bst, X.data[ind.val, , drop = FALSE])
    tmp <- mean(round5(pred2) != X.label[ind.val])
    if (tmp < (0.9 * error.save[i])) {
      error[i] <- tmp
      # impute with new (better) prediction
      pred <- predict(bst, X.data[indNA, , drop = FALSE])
      X2[indNA, i] <- round4(pred)
    }
  }
}
plot(error.save, error)
print(mean(error/ error.save, na.rm = TRUE))
error.save <- error

mean(round(X2[indNA.part]) != store.NA[seq_along(indNA.part)])
# 6.9 -> 7

#### ####
### first iteration:
# estimated:

# true:
# max_depth = 4, nrounds = 30, size = 50 -> 6.8 -> 6.74 -> 6.9 -> 6.8 -> 6.8
# with first component of SVD: 7%
# with gamma = 1 -> 6.66
# max_depth = 2 ->
# max_depth = 6, nrounds = 100, size = 50 ->
# max_depth = 4, nrounds = 30, size = 200 ->
mean(round(X2[indNA.part]) != store.NA[seq_along(indNA.part)])

plot(nbNA, error)
curve(10 / x, col = 2, lwd = 2, add = TRUE, from = 0)
curve(5 / x, col = 3, lwd = 2, add = TRUE, from = 0)
curve(2 / x, col = 4, lwd = 2, add = TRUE, from = 0)

popres.sub <- subset(popres, ind.col = 1:m.part)
snp_writeBed(popres.sub, "../../Téléchargements/plink_linux_x86_64/plink.bed")

# snp_readBed("../../Téléchargements/plink_linux_x86_64/popres_imputeBeagle.bed",
#             backingfile = "popresBeagle")
popres_beagle <- snp_attach("backingfiles/popresBeagle.rds")
X3 <- attach.BM(popres_beagle$genotypes)
mean(X3[indNA.part] != store.NA[seq_along(indNA.part)]) # 6.85
