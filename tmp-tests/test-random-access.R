x <- runif(1e9)
system.time(print(test1(x)))
# 2 sec
ind <- sample(seq_along(x) - 1L)
system.time(print(test2(x, ind)))
# 92 sec

microbenchmark::microbenchmark(
  test1(x),
  test2(x, ind),
  times = 10
)
