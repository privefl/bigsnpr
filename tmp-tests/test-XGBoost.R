library(xgboost)

data(agaricus.train, package='xgboost')
data(agaricus.test, package='xgboost')
train <- agaricus.train
test <- agaricus.test

bstSparse <- xgboost(data = train$data, label = train$label,
                     max_depth = 2, eta = 1, nthread = 2, nrounds = 2,
                     objective = "binary:logistic")

dtrain <- xgb.DMatrix(data = train$data, label = train$label)
bstDMatrix <- xgboost(data = dtrain, max_depth = 2, eta = 1,
                      nthread = 2, nrounds = 2,
                      objective = "binary:logistic",
                      verbose = 2)

pred <- predict(bstDMatrix, test$data)

# size of the prediction vector
print(length(pred))

# limit display of predictions to the first 10
print(head(pred))

require(bigsnpr)
print(AucSampleConf(pred, y = (test$label - 0.5) * 2))

prediction <- as.numeric(pred > 0.5)
print(head(prediction))
err <- mean(as.numeric(pred > 0.5) != test$label)
print(paste("test-error=", err)) # 0.02

require(LiblineaR)
svm <- LiblineaR(as.matrix(train$data), train$label)
pred2 <- predict(svm, as.matrix(test$data))
err2 <- mean(pred2$predictions != test$label)
print(paste("test-error=", err2)) # 0

dtrain <- xgb.DMatrix(data = train$data, label=train$label)
dtest <- xgb.DMatrix(data = test$data, label=test$label)
watchlist <- list(train=dtrain, test=dtest)

bst <- xgb.train(data=dtrain, max_depth=2, eta=1, nthread = 2,
                 nrounds=2, watchlist=watchlist,
                 objective = "binary:logistic")

bst <- xgb.train(data=dtrain, max_depth=2, eta=1, nthread = 2, nrounds=2,
                 watchlist=watchlist, eval_metric = "error",
                 eval_metric = "logloss", objective = "binary:logistic")

bst <- xgb.train(data=dtrain, booster = "gblinear", max_depth=2, nthread = 2, nrounds=2,
                 watchlist=watchlist, eval_metric = "error",
                 eval_metric = "logloss", objective = "binary:logistic")


importance_matrix <- xgb.importance(model = bst)
print(importance_matrix)
xgb.plot.importance(importance_matrix = importance_matrix)
plot(importance_matrix$Importance, type = "h")
