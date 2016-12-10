source("R/utils.R")
Rcpp::sourceCpp('src/pruning.cpp')

lims <- c(1, ncol(X))
S <- maf$S
exclude <- sample(ncol(X), size = 50)
size <- 500
X.chr <- X
ind.train <- seq(nrow(X.chr))
thr.corr <- 0.01

ind.chr <- seq2(lims)
m.chr <- length(ind.chr)
S.chr <- S[ind.chr]
ord.chr <- order(S.chr, decreasing = TRUE)
size <- min(size, m.chr)
s <- setdiff(-size:size, 0)

# init
keep <- rep(FALSE, m.chr)
remain <- rep(TRUE, m.chr)
remain[exclude] <- FALSE

for (ind in ord.chr) {
  if (remain[ind]) { # already excluded?
    ind.col <- intersect(ind + s, which(remain))

    res <- R_squared_chr(pBigMat = X.chr@address,
                         rowInd = ind.train,
                         colInd = ind.col,
                         colMat0 = X.chr[, ind])

    remain[ind.col[res > thr.corr]] <- FALSE
    keep[ind] <- TRUE
  }
}
ind.chr[keep]

