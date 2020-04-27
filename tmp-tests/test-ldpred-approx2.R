coef1 <- sapply(all_ldpred, function(.) .$beta_est * sqrt(N) * gwas$std.err)
pred.coef1 <- big_prodMat(G, coef1, ind.row = ind.val)
cor_pred <- cor(pred.coef1)
cor_pred[is.na(cor_pred)] <- 0
colMeans(cor_pred ** 2, na.rm = TRUE)

coef2 <- sapply(seq_along(all_ldpred), function(k) {

  h2 <- all_ldpred[[k]]$h2_est
  p <- all_ldpred[[k]]$p_est

  coeff <- N * h2 / m
  C <- coeff / (coeff + 1)
  C2 <- (coeff + 1) / (coeff + p)
  Vnc <- C / N * (1 - C * h2)
  Vc <- Vnc + C * h2 / (m * p)

  C2 * new_beta * sqrt(N) * gwas$std.err /
    (1 + (1 - p) / p * sqrt(Vc / Vnc) * exp(0.5 * new_beta ** 2 * (1 / Vc - 1 / Vnc)))
})
cor(coef1, coef2)
pred.coef2 <- big_prodMat(G, coef2, ind.row = ind.val)
diag(cor(pred.coef1, pred.coef2))
plot(pred.coef1[, 6], pred.coef2[, 6])

apply(pred.coef1, 2, function(pred)
  round(100 * drop(cor(pred, y2[ind.val])**2), 2))
apply(pred.coef2, 2, function(pred)
  round(100 * drop(cor(pred, y2[ind.val])**2), 2))
round(100 * drop(cor(pred, y2[ind.val])**2), 2)

plot(apply(pred.coef1, 2, function(pred)
  round(100 * drop(cor(pred, y2[ind.val])**2), 2)),
  colMeans(cor_pred ** 2))
