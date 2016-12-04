test <- snp_readHapMap3()

# ind.LD <- excludeLDreg(test)
# ind.keep <- snp_pruning(test, exclude = ind.LD)

X <- test$genotypes[,]
ind.NA <- which(is.na(X), arr.ind = TRUE)

ind <- seq(nrow(X))
Chr1 <- X[, 1:3000]
storage.mode(Chr1) <- "double"

ind.train <- sort(sample(ind, 600))
ind.test <- setdiff(ind, ind.train)

indNA <- sample(setdiff(1:3000, ind.NA[, 2]), 1)

require(xgboost)

dtrain <- xgb.DMatrix(data = Chr1[ind.train, -indNA],
                      label = Chr1[ind.train, indNA])
dtest <- xgb.DMatrix(data = Chr1[ind.test, -indNA],
                     label = Chr1[ind.test, indNA])
watchlist <- list(train=dtrain, test=dtest)

bst <- xgb.train(data = dtrain, max_depth=3, eta=0.3, nthread = 1,
                 nrounds=5, watchlist=watchlist, silent = 0,
                 objective = "multi:softprob", num_class = 3,
                 alpha = 10, lambda = 10)
importance_matrix <- xgb.importance(model = bst)
print(importance_matrix)
