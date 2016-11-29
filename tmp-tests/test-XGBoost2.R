require(bigsnpr)

popres <- AttachBigSNP("popres_sub1", "../thesis-celiac/popres/backingfiles/")
X <- popres$genotypes

ind <- which(popres$map$chromosome == 1)
Chr1 <- X[, ind]
storage.mode(Chr1) <- "double"

ind.train <- sort(sample(nrow(X), 1000))
ind.test <- setdiff(seq(nrow(X)), ind.train)

indNA <- 1000

dtrain <- xgb.DMatrix(data = Chr1[ind.train, -indNA],
                      label = Chr1[ind.train, indNA])
dtest <- xgb.DMatrix(data = Chr1[ind.test, -indNA],
                     label = Chr1[ind.test, indNA])
watchlist <- list(train=dtrain, test=dtest)

bst <- xgb.train(data = dtrain, max_depth=3, eta=1, nthread = 1,
                 nrounds=2, watchlist=watchlist, silent = 0,
                 objective = "multi:softmax", num_class = 3,
                 alpha = 1, lambda = 1)
# 1.81 with multi:softprob

pred <- predict(bst, dtest)
true <- getinfo(dtest, "label")
table(pred, true)
#plot(pred, true)
mean(round(pred) != true)

importance_matrix <- xgb.importance(model = bst)
print(importance_matrix)
xgb.plot.importance(importance_matrix = importance_matrix)
plot(importance_matrix$Feature,
     importance_matrix$Importance, type = "h")
summary(lm(Chr1[, indNA] ~ Chr1[, 1005]))
summary(lm(Chr1[, indNA] ~ Chr1[, 991]))
# Beware that indices begin at 0!!
