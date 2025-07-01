################################################################################

#' @importFrom bigassertr printf message2 warning2 stop2
#' @importFrom bigassertr assert_int assert_pos assert_01 assert_nona
#' @importFrom bigassertr assert_lengths assert_sorted assert_args
#' @importFrom bigassertr assert_noexist assert_exist assert_dir assert_ext
#' @importFrom bigassertr assert_type assert_class assert_package
#' @importFrom bigassertr assert_df_with_names assert_not_null
#' @importFrom bigparallelr assert_cores
assert_noNA <- bigstatsr:::assert_noNA

assert_bed <- function(x) {
  if (!(inherits(x, "bed") || inherits(x, "bed_light")))
    stop2("'%s' is not of class 'bed' or 'bed_light'.", deparse(substitute(x)))
}

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
        obj.bed     = "bigsnpr:::assert_bed(obj.bed)",
        infos.chr   = "assert_not_null(infos.chr)",
        # infos.pos   = "assert_pos(infos.pos)",
        ncores      = "assert_cores(ncores)",
        ind.row     = "assert_not_null(ind.row);   assert_int(ind.row);   assert_pos(ind.row)",
        ind.train   = "assert_not_null(ind.train); assert_int(ind.train); assert_pos(ind.train)",
        ind.col     = "assert_not_null(ind.col);   assert_int(ind.col);   assert_pos(ind.col)",
        ind.keep    = "assert_not_null(ind.keep);  assert_int(ind.keep);  assert_pos(ind.keep)",
        exclude     = "assert_int(exclude);   assert_pos(exclude)",
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
