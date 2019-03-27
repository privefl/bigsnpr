library(foreach)
N <- 1e4; K <- 200; L <- 50
system.time({
  scores_all_chr <- matrix(0, N, K)
  replicate(L, {
    all_prs_chr <- foreach(ic = 1:K, .combine = "cbind") %do% {
      rep(ic, N)
    }
    scores_all_chr <- scores_all_chr + all_prs_chr
  })
})

system.time({
  scores_all_chr2 <- replicate(K, rep(0, N), simplify = FALSE)
  replicate(L, {
    all_prs_chr2 <- foreach(ic = 1:K) %do% {
      rep(ic, N)
    }
    for (j in 1:K) {
      scores_all_chr2[[j]] <- scores_all_chr2[[j]] + all_prs_chr2[[j]]
    }
  })
  scores_all_chr2 <- do.call("cbind", scores_all_chr2)
})

all.equal(scores_all_chr, scores_all_chr2)
