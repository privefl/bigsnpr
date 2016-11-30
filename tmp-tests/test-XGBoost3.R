require(bigsnpr)

celiac <- AttachBigSNP("celiac_impute1_sub1")
X <- celiac$genotypes

ind <- which(celiac$fam$pheno == -1)
Chr1 <- X[, 1:3000]
storage.mode(Chr1) <- "double"

ind.train <- sort(sample(ind, 8000))
ind.test <- setdiff(ind, ind.train)

indNA <- runif(1, 1200, 1800)

dtrain <- xgb.DMatrix(data = Chr1[ind.train, -indNA],
                      label = Chr1[ind.train, indNA])
dtest <- xgb.DMatrix(data = Chr1[ind.test, -indNA],
                     label = Chr1[ind.test, indNA])
watchlist <- list(train=dtrain, test=dtest)

bst <- xgb.train(data = dtrain, max_depth=3, eta=0.3, nthread = 1,
                 nrounds=5, watchlist=watchlist, silent = 0,
                 objective = "multi:softprob", num_class = 3,
                 alpha = 1, lambda = 1)
# 14% with multi:softprob
# NEED TO ORDER -> define cusomized obj?

matrix(predict(bst, dtest)[1:18], 3)




pred <- predict(bst, dtest)
true <- getinfo(dtest, "label")
table(pred, true)
#plot(pred, true)
ind2 <- which(round(pred) != true)
print(ind3 <- ind.test[ind2])
celiac$fam$pop[ind3]
mean(round(pred) != true)


importance_matrix <- xgb.importance(model = bst)
print(importance_matrix)
xgb.plot.importance(importance_matrix = importance_matrix)
plot(importance_matrix$Feature,
     importance_matrix$Importance, type = "h", xlim = indNA + c(-10, 10))
summary(lm(Chr1[, indNA] ~ Chr1[, 9]))
summary(lm(Chr1[, indNA] ~ Chr1[, 991]))
# Beware that indices begin at 0!!
