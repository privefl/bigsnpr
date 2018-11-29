library(future)
library(foreach)

ncores <- 3
cl <- parallel::makeCluster(3)
avail <- bigstatsr::FBM(ncores, 1, type = "integer", init = 1)
doFuture::registerDoFuture()

res <- vector("list", 5)
for (i in seq_along(res)) {

  while (sum(avail[]) == 0) {
    cat("Waiting..\n")
    Sys.sleep(0.5)
  }
  ind.avail <- which(avail[] == 1)
  cat("Available:", length(ind.avail), "\n")

  plan(cluster, workers = cl[ind.avail])
  foo <- foreach(i = 3:1) %dopar% {
    Sys.sleep(i)
  }

  print(one <- ind.avail[1])
  avail[one] <- FALSE; print(avail[])
  res[[i]] <- cluster(workers = cl[one], {
    Sys.sleep(5)
    avail[one] <- TRUE
    i
  })
}

sapply(res, resolved)
parallel::stopCluster(cl)

