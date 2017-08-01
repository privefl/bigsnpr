code <- test <- readLines("https://raw.githubusercontent.com/privefl/bigstatsr/master/src/colstats.cpp")

fun_beg <- grep("template <class C>", code, fixed = TRUE)
fun_ends <- grep("^}$", code)
fun_end <- fun_ends[which(fun_ends > fun_beg)[1]]
code[fun_beg:fun_end]
tmp <- paste(code[fun_beg:fun_end], collapse = "\n")

source
