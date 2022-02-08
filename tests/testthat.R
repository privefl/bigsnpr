library(testthat)
library(bigsnpr)

for (k in 1:9)
  test_check("bigsnpr", filter = paste0(k, '-'))

# test_check("bigsnpr", filter = "BGEN")
