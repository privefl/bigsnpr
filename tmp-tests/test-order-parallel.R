library(foreach)

cl <- parallel::makeCluster(10)
doParallel::registerDoParallel(cl)


IND <- runif(20)
ORD <- order(IND, decreasing = TRUE)
ORD2 <- match(seq_along(ORD), ORD)


print(system.time(
  test <- foreach(ic = IND) %dopar% {
    Sys.sleep(ic)
    ic
  }
))

print(system.time(
  test2 <- foreach(ic = IND[ORD]) %dopar% {
    Sys.sleep(ic)
    ic
  }
))

foreach(ic = test2[ORD2], .combine = 'c') %do% {
  ic
}
unlist(test)


parallel::stopCluster(cl)
