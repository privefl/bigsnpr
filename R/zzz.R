.onLoad <- function(libname, pkgname) {
  options(bigsnpr.nploidy = 2)
}

.onUnload <- function(libpath) {
  options(bigsnpr.nploidy = NULL)
}
