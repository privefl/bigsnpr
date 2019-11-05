
x <- rep(c("01_88169_C_T", "1_88169_C_T"), 10000)
system.time(
  y <- gsubfn::strapply(
    X = x,
    pattern = "^(.+?)(_.+_.+_.+)$",
    FUN = function(x, y) paste0(ifelse(nchar(x) == 1, paste0("0", x), x), y),
    # perl = TRUE,
    # engine = "tcl",
    empty = stop("Wrong format of some SNPs."),
    simplify = 'c'
  )
)
# default: 15 sec
# perl:    25 sec
# tcl:     14 sec
# options(gsubfn.engine = "R"):

x2 <- rep(x, 100)
system.time(
  y2 <- sapply(strsplit(x, split = "_", fixed = TRUE), function(split) {
    if (length(split) != 4) stop("Wrong format of some SNPs.")
    if (nchar(split[1]) == 1) split[1] <- paste0("0", split[1])
    paste(split, collapse = "_")
  })
) # 0.2 sec
all.equal(y2, y)

system.time(
  y3 <- ifelse(grepl("^[[:digit:]]{1}_.*", x), paste0("0", x), x)
) # 0.02 sec
all.equal(y3, y)

# x <- rep(x, 100)
system.time(
  y4 <- ifelse(substr(x, 2, 2) == "_", paste0("0", x), x)
) # 0.1 sec
all.equal(y4, y)

system.time({
  y5 <- ifelse(substr(x, 2, 2) == "_", paste0("0", x), x)
  if (any(substr(y5, 3, 3) != "_")) stop("Wrong format of some SNPs.")
}) # 0.2 sec
all.equal(y5, y)
