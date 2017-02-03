require(bigsnpr)
require(xgboost)

celiac <- snp_attach("../thesis-celiac/backingfiles/celiac.bk")
X <- celiac$genotypes
X1 <- deepcopy(X, cols = which(celiac$map$chromosome == 1))


m <- ncol(X1)
breaks <- c(0, trunc(nrow(X) / c(1000, 200, 100)))
sizes <- c(10, 20, 30, 50)
nbNA <- integer(m)
error <- numeric(m)
all_ind <- seq(nrow(X))
ind.train <- sort(sample(all_ind, 8000))
ind.val <- setdiff(all_ind, ind.train)

interval <- function(i, size) {
  ind <- i + -size:size
  ind[ind >= 1 & ind <= m & ind != i]
}

round2 <- function(pred) {
  (pred > 0.5) + (pred > 1.5)
}

begin <- proc.time()
for (i in 1:m) {
  if (!(i %% 10)) print(i)
  X.label <- X[, i] * 1
  nbNA[i] <- l <- length(indNA <- which(is.na(X.label)))
  if (l > 0) {
    s <- sizes[max(which(l > breaks))]
    ind <- setdiff(ind.train, indNA)
    X.data <- X[, interval(i, s), drop = FALSE] * 1

    bst <- xgboost(data = X.data[ind, , drop = FALSE],
                   label = X.label[ind],
                   params = list(max_depth=5),
                   nrounds = 20,
                   nthread = 1, verbose = 0)
    pred <- predict(bst, X.data[indNA, , drop = FALSE])
    X1[indNA, i] <- round2(pred)
    pred2 <- predict(bst, X.data[ind.val, , drop = FALSE])
    print(error[i] <- mean(round2(pred2) != X.label[ind.val], na.rm = TRUE))
  }
}
print(proc.time() - begin)
# 208 sec for l > 0
# 30 sec with tradoff

plot(nbNA, error, cex = 0.5, pch = 19)
curve(10 / x, col = "red", add = TRUE,
      from = 1e-6, to = max(nbNA), n = 1e4)
curve(5 / x, col = "blue", add = TRUE,
      from = 1e-6, to = max(nbNA), n = 1e4)
curve(2 / x, col = "limegreen", add = TRUE,
      from = 1e-6, to = max(nbNA), n = 1e4)

nbNA2 <- rle(sort(nbNA, decreasing = TRUE))
nbNA3 <- rle(sort(nbNA))

imputeXGBoost <- function(X, nround = 20, max_depth = 5) {

}
