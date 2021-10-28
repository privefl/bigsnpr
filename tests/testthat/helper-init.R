################################################################################

library(testthat)
library(bigsnpr)

################################################################################

options(bigstatsr.check.parallel.blas = FALSE)

################################################################################

# https://github.com/hadley/testthat/issues/567
Sys.unsetenv("R_TESTS")

is_cran <- !identical(Sys.getenv("BIGSNPR_CRAN"), "false")
NCORES <- `if`(!is_cran && (parallel::detectCores() > 2) &&
                 identical(Sys.getenv("_R_CHECK_LIMIT_CORES_"), "false"), 2, 1)

################################################################################
