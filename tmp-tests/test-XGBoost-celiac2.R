celiac <- AttachBigSNP("celiac_impute1_sub1")
X <- celiac$genotypes
y <- celiac$fam$affection
pop <- celiac$fam$pop

MAX3 <- snp_MAX3(celiac)
ind.keep <- order(MAX3$pS)[1:1e4]

require(xgboost)


ind.train <- which(pop == 4)
ind.test <- which(pop == 2)

X.train <- X[ind.train, ind.keep] # 305 Mb
storage.mode(X.train) <- "double" # 610 Mb

dtrain <- xgb.DMatrix(data = X.train, label = y[ind.train] - 1)
rm(X.train)

X.test <- X[ind.test, ind.keep] # 278 Mb
storage.mode(X.test) <- "double" # 557 Mb

dtest <- xgb.DMatrix(data = X.test, label = y[ind.test] - 1)
rm(X.test)

watchlist <- list(train = dtrain, test = dtest)

bst <- xgb.train(params = list(max_depth = 10, colsample_bylevel = 0.8,
                               objective = "binary:logitraw",
                               alpha = 5, lambda = 100),
                 data = dtrain, watchlist = watchlist,
                 nthread = 3, nrounds = 100, booster = "dart")

importance_matrix <- xgb.importance(model = bst)
print(importance_matrix)
print(ind <- match(1:100, importance_matrix$Feature))
plot(importance_matrix$Feature[1:100],
     importance_matrix$Importance[1:100], type = "h")
plot(MAX3$pS[ind.keep[as.numeric(importance_matrix$Feature)]],
     importance_matrix$Importance, pch = 19, cex = 0.5, log = "xy")



K <- 20
ind2 <- as.numeric(importance_matrix$Feature[1:K])
X.train <- X[ind.train, ind.keep[ind2]] # 305 Mb
storage.mode(X.train) <- "double" # 610 Mb

dtrain2 <- xgb.DMatrix(data = X.train, label = y[ind.train] - 1)
rm(X.train)

X.test <- X[ind.test, ind.keep[ind2]] # 278 Mb
storage.mode(X.test) <- "double" # 557 Mb

dtest2 <- xgb.DMatrix(data = X.test, label = y[ind.test] - 1)
rm(X.test)

watchlist2 <- list(train = dtrain2, test = dtest2)

bst2 <- xgb.train(params = list(max_depth = 10, colsample_bylevel = 1,
                                objective = "binary:logitraw",
                                alpha = 5, lambda = 100),
                  data = dtrain2, watchlist = watchlist2,
                  nthread = 3, nrounds = 100)
