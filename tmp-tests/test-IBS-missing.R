require(bigsnpr)
celiac <- snp_attach("backingfiles/celiac300_sub1.rds")
G <- celiac$genotypes
G@code <- c(0:2, NA, 0:2, seq(0, 2, by = 0.01), rep(NA, 48))

# tmpfile <- paste0("tmp_transpose", sub("0\\.", "", runif(1)))
# Gt <- big_transpose(G, backingfile = tmpfile,
#                     backingpath = "backingfiles")


print(system.time(
  test <- big_tcrossprodSelf(G, fun.scaling = snp_scaleBinom())
)) # 2 min for 15283 x 16325

test$K[1:5, 1:5]


print(system.time(
  test2 <- big_randomSVD(G, fun.scaling = snp_scaleBinom(),
                         k = 50, ncores = 4, verbose = TRUE,
                         ind.col = seq(1, ncol(G), by = 100))
)) # 4 min for 15283 x 16325
# need pruning absolutely

U <- test2$u %*% diag(test2$d)
norms <- sqrt(rowMeans(U^2))
K <- tcrossprod(U)
Rcpp::sourceCpp('tmp-save/center.cpp')
K <- toCorr(K, norms)
K <- K / 50
K[1:5, 1:5]

fun_IBS <- function(X.,
                    ind.row = rows_along(X.),
                    ind.col = cols_along(X.)) {
  m <- length(ind.col)
  data.frame(mean = rep(1, m), sd = rep(1, m))
} # code - 1 would do the same thing

print(system.time(
  test32 <- big_tcrossprodSelf(G, fun.scaling = fun_IBS)
)) # 2 min for 15283 x 16325

test32$K[1:5, 1:5]


d <- data.table::fread("../../Téléchargements/plink_linux_x86_64/plink.genome",
                       data.table = FALSE)

require(dplyr)
d$ind1 <- match(d$IID1, celiac$fam$sample.ID)
d$ind2 <- match(d$IID2, celiac$fam$sample.ID)
d$cos <- test32$K[as.matrix(d[, c("ind1", "ind2")])]
plot2 <- function(x, y) {
  plot(x, y)
  abline(lm(y ~ x), col = "red")
}
plot2(d$PI_HAT, d$cos)

print(system.time(
  test.dist <- dist(test2$u)
))


library(ggplot2)
d2 <- subset(d, PI_HAT > 0.1)
qplot(Z0,Z1,data=d2,colour=RT, cex = PI_HAT)

subset(d2, PI_HAT > 0.185)
ind <- match(c("74231_H09_BLOOD291884", "WG0012843-DNA_B04_7014"), celiac$fam$sample.ID)
twoIndiv <- t(attach.BM(G)[ind, ])



miss <- data.table::fread("../../Téléchargements/plink_linux_x86_64/plink.missing",
                       data.table = FALSE)

plot(miss$F_MISS_A, miss$F_MISS_U, cex = 0.6, pch = 19,
     col = (miss$P < 1e-8) + (miss$P < 1e-4) + (miss$P < 0.01) + (miss$P < 0.05) + 1); abline(0, 1, col = "red", lwd = 2)


N.ca <- 4533
N.co <- 10750

nbNA.ca <- round(miss$F_MISS_A * N.ca)
nbNA.co <- round(miss$F_MISS_U * N.co)

arr <- array(dim = c(2, 2, length(nbNA.ca)))
arr[1, 1, ] <- nbNA.ca
arr[2, 1, ] <- nbNA.co
arr[1, 2, ] <- N.ca - nbNA.ca
arr[2, 2, ] <- N.co - nbNA.co
print(system.time(
  p.fisher <- apply(arr, 3, function(mat) fisher.test(mat)$p.value)
)) # 90 sec

storage.mode(arr) <- "integer"
print(system.time(
  p.fisher2 <- apply(arr, 3, fisher)
)) # 8 sec -> 5 -> 4.5 -> 4
all.equal(p.fisher, p.fisher2)

indAuto <- which(celiac300$map$chromosome <= 22)
miss <- miss[indAuto, ]
plot(miss$F_MISS_A, miss$F_MISS_U, cex = 0.6, pch = 19,
     col = (miss$P < 1e-8) + (miss$P < 1e-4) + (miss$P < 0.01) + (miss$P < 0.05) + 1); abline(0, 1, col = "red", lwd = 2)
