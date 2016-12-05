require(bigsnpr)
require(xgboost)

celiac <- AttachBigSNP("celiac")
X <- celiac$genotypes
X.sub <- X[, 1:1000]
storage.mode(X.sub) <- "double"
breaks <- c(0, trunc(nrow(X) / c(1000, 200, 100)))
sizes <- c(10, 20, 30, 50)
indNA <- which(is.na(X.sub), arr.ind = TRUE)
listNA <- apply(is.na(X.sub), 2, which)
nbNA <- integer(1000)
error <- numeric(1000)
params <- list(max_depth=5, colsample_bylevel = 1)
all_ind <- seq(nrow(X))
ind.train <- sort(sample(all_ind, 8000))
ind.val <- setdiff(all_ind, ind.train)

interval <- function(size) setdiff(-size:size, 0)

begin <- proc.time()
for (i in 51:150) {
  nbNA[i] <- l <- length(indNA <- listNA[[i]])
  if (l > 0) {
    s <- interval(sizes[max(which(l > breaks))])
    ind <- setdiff(ind.train, indNA)
    bst <- xgboost(data = X.sub[ind, i + s], label = X.sub[ind, i],
                   params = params, nrounds = 20, nthread = 1, verbose = 0)
    pred <- predict(bst, X.sub[indNA, i + s, drop = FALSE])
    pred2 <- predict(bst, X.sub[ind.val, i + s, drop = FALSE])
    print(error[i] <- mean(round(pred2) != X.sub[ind.val, i], na.rm = TRUE))
  }
}
print(proc.time() - begin)
# 208 sec for l > 0
# 30 sec with tradoff

nbNA2 <- rle(sort(nbNA, decreasing = TRUE))
nbNA3 <- rle(sort(nbNA))


assess3 <- function(indNA, param, ...) {
  bst <- xgboost(params = param, data = Chr1[ind.train, indNA + s2],
                 label = Chr1[ind.train, indNA],
                 nthread = 1, ...)
  pred <- predict(bst, Chr1[ind.test, indNA + s2])
  true <- Chr1[ind.test, indNA]
  mean(round(pred) != true)
}

require(foreach)
res <- foreach(i = 1:20, .combine = 'cbind') %do% {
  indNA <- round(runif(1, 1200, 1800))
  time1 <- system.time(
    xgb1 <- assess3(indNA, ,
                    )
