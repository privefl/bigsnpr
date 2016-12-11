celiac <- AttachBigSNP("celiac_impute1_sub1")
X <- celiac$genotypes
y <- celiac$fam$affection

require(xgboost)

ind.train <- sort(sample(nrow(X), 1e4))
ind.test <- setdiff(1:nrow(X), ind.train)

MAX3 <- snp_MAX3(celiac, ind.train = ind.train)
ind.keep <- order(MAX3$pS)[1:1e4]

X.train <- X[ind.train, ind.keep] # 305 Mb
storage.mode(X.train) <- "double" # 610 Mb

dtrain <- xgb.DMatrix(data = X.train, label = y[ind.train] - 1)
rm(X.train)

X.test <- X[ind.test, ind.keep] # 278 Mb
storage.mode(X.test) <- "double" # 557 Mb

dtest <- xgb.DMatrix(data = X.test, label = y[ind.test] - 1)
rm(X.test)

watchlist <- list(train = dtrain, test = dtest)


ratio <- sum(y == 1) / sum(y == 2)
bst <- xgb.train(params = list(max_depth = 6, colsample_bylevel = 1,
                               objective = "binary:logitraw",
                               alpha = 10, lambda = 100, gamma = 1,
                               scale_pos_weight = ratio),
                 data = dtrain, watchlist = watchlist,
                 nthread = 3, nrounds = 100,
                 early_stopping_rounds = 5)
xgb.plot.tree(model = bst, n_first_tree = 0)

importance_matrix <- xgb.importance(model = bst)
print(importance_matrix)
print(ind <- match(1:100, importance_matrix$Feature))
# plot(MAX3$S[ind.keep[as.numeric(importance_matrix$Feature)]],
#      importance_matrix$Gain, pch = 19, cex = 0.5, log = "xy")

plot(predict(bst, dtest), jitter(y[ind.test]),
     pch = 10, cex = 0.5, col = pop[ind.test])
plot(predict(bst, dtrain), jitter(y[ind.train]),
     pch = 10, cex = 0.5, col = pop[ind.train])


model <- xgb.dump(bst, with_stats = TRUE)
model[1:10]
xgb.plot.tree(model = bst, n_first_tree = 0)

for (i in 1:4) {
  ind <- as.numeric(importance_matrix$Feature[i])
  print(table(X[, ind.keep[ind]], y))
}

for (i in 1:10) {
  print(table(X[, ind.keep[i]], y))
}
