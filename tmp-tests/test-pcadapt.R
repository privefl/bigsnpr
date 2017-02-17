#' Date: 2017-02-17
#' Object: Test pcadapt with pruning or not
#' Results: Removing LRLDR doesn't make them disappear in scans.


require(bigstatsr)
require(bigsnpr)
require(foreach)
require(robust)

popres <- snp_attach("../thesis-celiac/popres/backingfiles/popres.bk")
X <- popres$genotypes

fregion<-function(x)
{
  if((x=="Portugal" )|(x=="Spain"))
    return ("SW Europe")
  else if((x=="Greece") | (x=="Turkey") | (x=="Cyprus") | (x=="Macedonia,"))
    return("SE Europe")
  else if((x=="Italy"))
    return("Italy")
  else if((x=="?orway") |  (x=="Finland") | (x=="Sweden") | (x=="Denmark"))
    return("FennoScandia")
  else if((x=="France") |  (x=="Swiss-French") | (x=="Swiss-German") | (x=="Swiss-Italian") | (x=="Belgium"))
    return("Western Europe")
  else if(  (x=="Czech")  |(x=="Romania")|(x=="Bulgaria")| (x=="Hungary")| (x=="Albania") | (x=="Macedonia")| (x=="Slovakia"))
    return("Eastern Europe")
  else if ((x=="Croatia") | (x=="Serbia")|(x=="Bosnia")| (x=="Slovenia")| (x=="Kosovo") )
    return("Former Yugoslavia")
  else if((x=="United")  |  (x=="Scotland")  |  (x=="Ireland"))
    return("Anglo-Irish Isles")
  else if((x=="Germany")  |  (x=="Austria") |  (x=="Poland") | (x=="?etherlands"))
    return("Central Europe")
  else if((x=="Russian")  |  (x=="Ukraine") |  (x=="Latvia") )
    return ("Former USSR")
  else return ("problem")
}
popres$fam$pop = sapply(popres$fam$family.ID, fregion, USE.NAMES = FALSE)

n <- nrow(X)
m <- ncol(X)
maf <- snp_MAF(X)
maf.NOK <- (maf < 0.05)
in.LRLDR <- (1:m %in% snp_indLRLDR(popres))

pruning <- 1

if (pruning == 2) { ### With stringent pruning
  ind.keep <- snp_clumping(popres, S = maf, thr.corr = 0.05,
                           exclude = which(maf.NOK))
  X.svd <- big_randomSVD(X, fun.scaling = snp_scaleBinom(),
                         ind.col = ind.keep)
} else if (pruning == 0) { ### Without pruning
  X.svd <- big_randomSVD(X, fun.scaling = snp_scaleBinom(),
                         ind.col = which(!maf.NOK))
} else if (pruning == 1) { ### Without the long-range LD regions
  ind.keep <- snp_clumping(popres, S = maf, thr.corr = 0.5,
                           exclude = which(maf.NOK | in.LRLDR))
  X.svd <- big_randomSVD(X, fun.scaling = snp_scaleBinom(),
                         ind.col = ind.keep)
} else {
  stop("You shouldn't be here.")
}


plot(X.svd$u[, 3:4], col = as.factor(popres$fam$pop), pch = 19, cex = 0.5)

ms <- big_colstats(X)
z.scores <- linRegPcadapt(X@address, U = X.svd$u[, 1:5], rowInd = seq(n),
                          center = ms$sum / n, scale = sqrt(ms$var))
# z.scores2 <- foreach(ic = 1:3, .combine = 'cbind') %do% {
#   big_univRegLin(X, X.svd$u[, ic])$t.score
# } # equivalent

distCovRob <- rep(NA_real_, m)
are.NA <- (rowSums(is.na(z.scores)) > 0 | maf.NOK) # all MAF < 5% -> NA when MAF ~ 0
distCovRob[!are.NA] <- covRob(z.scores[!are.NA, ], estim = "pairwiseGK")$dist

cols <- ifelse(in.LRLDR, 3, (popres$map$chromosome %% 2) + 1)
plot(distCovRob, pch = 19, cex = 0.5, col = cols)

ind <- which(distCovRob > 50)
print(ind %in% ind.keep)

test <- big_apply(X, FUN = function(x, ind, U, n) {
  ms <- big_colstats(x, ind.col = ind)
  x.part <- sub.big.matrix(x, firstCol = head(ind, 1), lastCol = tail(ind, 1))
  z.scores <- linRegPcadapt(x.part@address, U, rowInd = seq(n),
                            center = ms$sum / n, scale = sqrt(ms$var))
}, .combine = 'rbind', block.size = m, ind.arg = TRUE, ncores = 6,
U = X.svd$u[, 1:5], n = n)
print(all.equal(test, z.scores))

distCovRob2 <- subset(distCovRob, !is.na(distCovRob))
print(lamGC <- median(distCovRob2) / qchisq(0.5, df = 5, lower.tail = FALSE))
pS <- pchisq(distCovRob2 / lamGC, df = 5, lower.tail = FALSE)
hist(pS)
m2 <- length(pS)
plot(-log10(m2:1 / (m2 + 1)), sort(-log10(pS)), pch = 19, cex = 0.5)
abline(0, 1, col = "red")

plot(-log10(pchisq(distCovRob / lamGC, df = 5, lower.tail = FALSE)),
     pch = 19, cex = 0.5, col = cols, ylab = "-log10(p-value)")
