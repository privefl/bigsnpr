require(bigsnpr)

celiac <- AttachBigSNP("celiac_impute1_sub1")
X <- celiac$genotypes

ind <- which(celiac$fam$pheno == -1)
Chr1 <- X[, 1:3000]
storage.mode(Chr1) <- "double"

ind.train <- sort(sample(ind, 8000))
ind.test <- setdiff(ind, ind.train)

size <- 200
s <- setdiff(-size:size, 0)

assess <- function(indNA, param, ...) {
  dtrain <- xgb.DMatrix(data = Chr1[ind.train, indNA + s],
                        label = (Chr1[ind.train, indNA]))
  dtest <- xgb.DMatrix(data = Chr1[ind.test, indNA + s],
                       label = (Chr1[ind.test, indNA]))
  watchlist <- list(train=dtrain, test=dtest)

  bst <- xgb.train(param, data = dtrain, nthread = 1,
                   watchlist=watchlist, ...)

  print(indNA)
  importance_matrix <- xgb.importance(model = bst)
  print(importance_matrix)
  xgb.plot.importance(importance_matrix = importance_matrix)

  pred <- predict(bst, dtest)
  pred2 <- round(pred)
  true <- getinfo(dtest, "label")
  table(pred2, true)
  mean(round(pred) != true)
}

assess2 <- function(indNA, param, ...) {
  bst <- xgboost(params = param, data = Chr1[ind.train, indNA + s],
                 label = Chr1[ind.train, indNA],
                 nthread = 1, ...)
  pred <- predict(bst, Chr1[ind.test, indNA + s])
  true <- Chr1[ind.test, indNA]
  mean(round(pred) != true)
}

require(foreach)
res <- foreach(i = 1:20, .combine = 'cbind') %do% {
  indNA <- round(runif(1, 1200, 1800))
  time1 <- system.time(
    xgb1 <- assess2(indNA, param = list(max_depth=5, eta=0.3,
                                        alpha = 10, lambda = 10),
                    verbose = 0, nrounds = 100)
  )[3]
  time2 <- system.time(
    xgb2 <- assess2(indNA, param = list(max_depth=5, eta=0.3,
                                        alpha = 20, lambda = 20),
                    verbose = 0, nrounds = 100)
  )[3]
  c(xgb1, xgb2, time1, time2)
}
plot(t(res))
abline(0, 1, col = "red")
rowMeans(res)
# reg:linear faster and maybe better than multi:softmax
# booster = "gblinear" is faster but worst than the default gbtree
# alpha = 1000 is baaaad
# alpha = 10 is better
# max_depth = 5 is better
# eta = 0.3 is better than 1
