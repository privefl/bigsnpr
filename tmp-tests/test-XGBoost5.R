
dtrain.mat <- Chr1[ind.train, indNA + s]
m <- ncol(dtrain.mat)
preds <- c(1/3, 2/3, rep(0, m))

obj <- function(preds, dtrain) {
  labels <- getinfo(dtrain, "label")
  ind1 <- which(labels == 0)
  ind2 <- which(labels == 1)
  ind3 <- which(labels == 2)

  a <- preds[1]
  b <- preds[2]
  w <- preds[-(1:2)]

  e_a <- exp(a)
  e_b <- exp(b)

  # y1
  mat1 <- dtrain.mat[ind1, ]
  e_xw1 <- exp(mat1 %*% w)
  e_axw1 <- e_a + e_xw1
  C1 <- e_xw1 / e_axw1
  C2 <- e_a * e_xw1 / e_axw1^2
  f1.w <- crossprod(mat1, C1)
  f1.w2 <- crossprod(mat1^2, C2)
  f1.a <- - sum(C1)
  f1.a2 <- sum(C2)

  # y2
  mat2 <- dtrain.mat[ind2, ]
  e_xw2 <- exp(mat2 %*% w)
  e_axw2 <- e_a + e_xw2
  e_bxw2 <- e_b + e_xw2
  f2.w <- crossprod(mat2, (e_xw2^2 - e_a * e_b) / (e_axw2 * e_bxw2))
  f2.w2 <- crossprod(mat2^2, e_xw2 / (e_axw2 * e_bxw2)^2 *
                       (4 * e_xw2 * e_a * e_b + (e_b + e_a) * e_xw2^2 +
                          e_a * e_b^2 + e_b * e_a^2))
  C1 <- (e_a - e_b) * e_axw2
  f2.a <- - sum(e_bxw2 * e_a / C1)
  f2.a2 <- sum(e_bxw2 * e_a * (e_a^2 + e_bxw2) / C1^2)
  C2 <- (e_b - e_a) * e_bxw2
  f2.b <- - sum(e_axw2 * e_b / C2)
  f2.b2 <- sum((e_axw2 * e_b * (e_b^2 + e_axw2)) / C2^2)

  # y3
  mat3 <- dtrain.mat[ind3, ]
  e_xw3 <- exp(mat3 %*% w)
  e_bxw3 <- e_b + e_xw3
  C1 <- e_b / e_bxw3
  C2 <- e_xw3 * e_b / e_bxw3^2
  f3.w <- - crossprod(mat3, C1)
  f3.w2 <- crossprod(mat3^2, C2)
  f3.b <- sum(C1)
  f3.b2 <- sum(C2)

  f.w <- f1.w + f2.w + f3.w
  f.w2 <- f1.w2 + f2.w2 + f3.w2
  f.a <- f1.a + f2.a
  f.a2 <- f1.a2 + f2.a2
  f.b <- f2.b + f3.b
  f.b2 <- f2.b2 + f3.b2

  list(grad = c(f.a, f.b, f.w), hess = c(f.a2, f.b2, f.w2))
}

test <- obj(preds, dtrain)
preds <- preds - test$grad / test$hess
