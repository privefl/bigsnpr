################################################################################

library(testthat)
library(bigsnpr)
library(bigparallelr)

################################################################################

options(bigstatsr.check.parallel.blas = FALSE)

################################################################################

# https://github.com/hadley/testthat/issues/567
Sys.unsetenv("R_TESTS")

not_cran <- identical(Sys.getenv("BIGSNPR_CRAN"), "false")
NCORES <- `if`(not_cran, 2, 1)

################################################################################
