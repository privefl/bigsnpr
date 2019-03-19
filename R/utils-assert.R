################################################################################

check_args <- function(...) {

  if (getOption("bigstatsr.check.args")) {

    args <- as.list(parent.frame())

    check <- c(
      list(...),  # possible to "overwrite" following defaults
      list(
        x           = "assert_class(x, 'bigSNP')",
        G           = "assert_class(G, 'FBM.code256'); assert_noNA(G)",
        Gna         = "assert_class(Gna, 'FBM.code256')",
        infos.chr   = "assert_int(infos.chr); assert_pos(infos.chr); assert_sorted(infos.chr)",
        infos.pos   = "assert_int(infos.pos); assert_pos(infos.pos); assert_sorted(infos.pos)",
        ncores      = "assert_cores(ncores)",
        ind.row     = "assert_int(ind.row);   assert_pos(ind.row)",
        ind.train   = "assert_int(ind.train); assert_pos(ind.train)",
        ind.col     = "assert_int(ind.col);   assert_pos(ind.col)",
        ind.keep    = "assert_int(ind.keep);  assert_pos(ind.keep)",
        fun.scaling = "assert_args(fun.scaling, c('ind.row', 'ind.col'))",
        gwas        = "assert_class(gwas, 'mhtest')",
        y01.train   = "assert_01(y01.train)"
      )
    )

    for (i in match(names(args), names(check)))
      if (!is.na(i)) with(args, eval(parse(text = check[i])))
  }
}

################################################################################

# MISSING VALUES
assert_noNA        <- bigstatsr:::assert_noNA
assert_nona        <- bigstatsr:::assert_nona
# DIRECTORY
assert_dir         <- bigstatsr:::assert_dir
# ARGS
assert_args        <- bigstatsr:::assert_args
# NUMBER OF CORES
assert_cores       <- bigstatsr:::assert_cores
# INTEGERS
assert_int         <- bigstatsr:::assert_int
# POSITIVE INDICES
assert_pos         <- bigstatsr:::assert_pos
# 0s AND 1s
assert_01          <- bigstatsr:::assert_01
# CLASS
assert_class       <- bigstatsr:::assert_class
# FILE EXISTS
assert_exist       <- bigstatsr:::assert_exist
assert_noexist     <- bigstatsr:::assert_noexist
# LENGTH
assert_lengths     <- bigstatsr:::assert_lengths

# EXTENSION
assert_ext <- function(file, ext) {
  if (!grepl(sprintf("\\.%s$", ext), file))
    stop2("Extension of '%s' must be '.%s'.", file, ext)
}

# TYPEOF
assert_type <- function(x, type)  {
  if (typeof(x) != type)
    stop2("'%s' is not of type '%s'.", deparse(substitute(x)), type)
}

################################################################################

# SORT
assert_sorted <- function(x)  {
  if (is.unsorted(x))
    stop2("'%s' is not sorted.", deparse(substitute(x)))
}

################################################################################
