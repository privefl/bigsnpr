celiac <- AttachBigSNP("celiac_impute1_sub1")

MAX3 <- snp_MAX3(celiac)

lam <- median(MAX3$S) / qchisq(0.5, df = 1)
S2 <- MAX3$S / lam

m <- length(S2)
expect <- 1:m / (m + 1)
lpS2 <- -log10(pchisq(S2, df = 1, lower.tail = FALSE))

plot(-log10(rev(expect)), sort(lpS2), cex = 0.5, pch = 19)
abline(0, 1, col = "red")

plot(lpS2, col = celiac$map$chromosome, pch = 19, cex = 0.5)


# without chr6
S3 <- MAX3$S
S3[celiac$map$chromosome == 6] <- NA
q <- 0.9
lam3 <- quantile(S3, probs = q, na.rm = TRUE) /
  qchisq(q, df = 1)
S4 <- S3 / lam3
lpS4 <- -log10(pchisq(S4, df = 1, lower.tail = FALSE))

m2 <- m - sum(is.na(lpS4))
expect2 <- 1:m2 / (m2 + 1)
plot(-log10(rev(expect2)), sort(lpS4), cex = 0.5, pch = 19)
abline(0, 1, col = "red")

plot(lpS4, col = celiac$map$chromosome, pch = 19, cex = 0.5)
