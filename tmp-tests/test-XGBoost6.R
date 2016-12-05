preds <- rnorm(ind.train)

obj <- function(preds, dtrain) {
  labels <- getinfo(dtrain, "label")
  ind1 <- which(labels == 0)
  ind2 <- which(labels == 1)
  ind3 <- which(labels == 2)

  n <- length(labels)
  grad <- numeric(n)
  hess <- numeric(n)

  m <- mean(labels)
  a <- 0.5 - m
  b <- 1.5 - m
  e_a <- exp(a)
  e_b <- exp(b)
  e_ab <- e_a * e_b

  # y1
  e_x1 <- exp(preds[ind1])
  C1 <- e_x1 + e_a
  grad[ind1] <- e_x1 / C1
  hess[ind1] <- e_a * e_x1 / C1^2

  # y2
  e_x2 <- exp(preds[ind2])
  C1 <- e_x2 + e_a
  C2 <- e_x2 + e_b
  grad[ind2] <- (e_x2^2 - e_ab) / (C1 * C2)
  hess[ind2] <- e_x2 * (4 * e_x2 * e_ab + (e_b + e_a) * (e_x2^2 + e_ab)) /
    (C1 * C2)^2

  # y3
  e_x3 <- exp(preds[ind3])
  C1 <- e_x3 + e_b
  grad[ind3] <- - e_b / C1
  hess[ind3] <- e_b * e_x3 / C1^2

  list(grad = grad, hess = hess)
}

phi <- function(theta, x) {
  1 / (1 + exp(x - theta))
}




feval <- function(preds, dtrain) {
  labels <- getinfo(dtrain, "label")
  pred1 <- preds[labels == 0]
  pred2 <- preds[labels == 1]
  pred3 <- preds[labels == 2]
  logLik <- sum(log(phi(a, pred1))) +
    sum(log(phi(b, pred2) - phi(a, pred2))) +
    sum(log(1 - phi(b, pred2)))
  list(metric = 'negative-logLik', value = -logLik)
}


feval2 <- function(preds, dtrain) {
  pred <- (preds > a) + (preds > b)
  list(metric = 'accuracy', value = mean(pred == labels))
}

test <- obj(preds, dtrain)
preds <- preds - 0.5 * test$grad / test$hess
plot(preds)

plot(labels, preds)


indNA <- round(runif(1, 1200, 1800))
ind.train <- sort(sample(nrow(X), 9000))
ind.remains <- setdiff(1:nrow(X), ind.train)
ind.val <- sort(sample(ind.remains, 3000))
ind.test <- setdiff(ind.remains, ind.val)
dtrain <- xgb.DMatrix(data = Chr1[ind.train, indNA + s],
                      label = (Chr1[ind.train, indNA]))
labels <- getinfo(dtrain, "label")
m <- mean(labels)
a <- 0.5 - m
b <- 1.5 - m
dtest <- xgb.DMatrix(data = Chr1[ind.test, indNA + s],
                     label = (Chr1[ind.test, indNA]))
dval <- xgb.DMatrix(data = Chr1[ind.val, indNA + s],
                     label = (Chr1[ind.val, indNA]))
watchlist <- list(train=dtrain, test=dtest)

bst2 <- xgb.train(params = list(max_depth=5, colsample_bylevel = 0.5),
                  data = dtrain, nthread = 1,
                  watchlist=watchlist, nrounds = 50)

bst <- xgb.train(params = list(max_depth=5, colsample_bylevel = 0.5),
                 data = dtrain, nthread = 1, obj = obj,
                 feval = feval, maximize = TRUE,
                 watchlist=watchlist, nrounds = 50)


pred <- predict(bst, dtest)
true <- getinfo(dtest, "label")
plot(predict(bst2, dtest), pred, col = true + 1, pch = 19, cex = 0.5)
abline(h = c(-m, 2 - m))
mean(round(predict(bst2, dtest)) != true)
pred2 <- predict(bst, dtrain)
true2 <- getinfo(dtrain, "label")
#points(predict(bst2, dtrain), pred2, col = true2 + 4, pch = 19, cex = 0.5)
abline(v = c(0.5, 1.5))
mean(round(predict(bst2, dtrain)) != true2)

true3 <- getinfo(dval, "label")
mean(round(predict(bst2, dval)) != true3)



plot(predict(bst2, dtest), pred, col = celiac$fam$pop[ind.test], pch = 19, cex = 0.5)


m <- mean(true2)
a <- 0.5 - m
b <- 1.5 - m
# print(mean(round(pred2) != true2))
# mean(round(pred) != true)

trans <- function(pred) {
  (pred > a) + (pred > b)
}
print(mean(trans(pred2) != true2))
mean(trans(pred) != true)
table(trans(pred), true)

plot(density(pred[true == 0]), xlim = c(-3,2), ylim = c(0, 0.6))
lines(density(pred[true == 1]), col = 2)
lines(density(pred[true == 2]), col = 3)

