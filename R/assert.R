# shorter name for function `assert_that`
assert <- assertthat::assert_that


# TYPE
has_type <- function(x, type) typeof(x) == type

assertthat::on_failure(has_type) <- function(call, env)
  sprintf("%s is not of type %s", deparse(call$x), deparse(call$type))
