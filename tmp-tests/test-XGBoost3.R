require(bigsnpr)

celiac <- AttachBigSNP("celiac_sub2_impute1",
                       backingpath = "../thesis-celiac/backingfiles/")
X <- celiac$genotypes

ind <- which(celiac$fam$pheno == -1)
Chr1 <- X[, 1:3000]
storage.mode(Chr1) <- "double"

ind.train <- sort(sample(ind, 7000))
ind.test <- setdiff(ind, ind.train)

size <- 200
s <- setdiff(-size:size, 0)

require(xgboost)
assess <- function(indNA, param, ...) {
  dtrain <- xgb.DMatrix(data = Chr1[ind.train, indNA + s],
                        label = (Chr1[ind.train, indNA]))
  dtest <- xgb.DMatrix(data = Chr1[ind.test, indNA + s],
                       label = (Chr1[ind.test, indNA]))
  watchlist <- list(train=dtrain, test=dtest)

  bst <- xgb.train(param, data = dtrain, nthread = 1,
                   watchlist=watchlist, ...)

  pred <- predict(bst, dtest)
  true <- getinfo(dtest, "label")
  pred2 <- predict(bst, dtrain)
  true2 <- getinfo(dtrain, "label")
  print(mean(round(pred2) != true2))
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
    xgb1 <- assess2(indNA, param = list(max_depth=5,
                                        colsample_bylevel = 1),
                    verbose = 0, nrounds = 50)
  )[3]
  time2 <- system.time(
    xgb2 <- assess2(indNA, param = list(max_depth=5, colsample_bylevel = 0.5),
                    verbose = 0, nrounds = 50)
  )[3]
  c(xgb1, xgb2, time1, time2)
}
plot(t(res))
abline(0, 1, col = "red")
rowMeans(res)
# reg:linear faster and maybe better than multi:softmax
# booster = "gblinear" is faster but worst than the default gbtree
# alpha = 1000 is baaaad
# alpha = 10 is better and even alpha = 0 -> 0.031 vs 0.286
# max_depth = 5 is better
# eta = 0.3 is better than 1
# 0.0364 vs 0.0357 for nrounds = 50 vs 200, but 2 times slower
# 0.051 vs 0.062 for nrounds = 50 vs 10, but 4 times slower
# min_child_weight is OK
# subsample = 1 is OK

# print(system.time(
#   bst <- xgboost(params = list(), data = Chr1[ind.train, indNA + s],
#                  label = Chr1[ind.train, indNA],
#                  nthread = 1, nrounds = 100)
# ))
#
# lo <- c(TRUE, FALSE)
# print(system.time(
#   bst <- xgboost(params = list(), data = Chr1[ind.train, indNA + s][lo, ],
#                  label = Chr1[ind.train, indNA][lo],
#                  nthread = 1, nrounds = 100)
# )) # algo seems to be linear with nb of indiv

print(system.time(print(
  assess2(indNA, param = list(colsample_bylevel = 0.5),
         verbose = 0, nrounds = 20)
)))


