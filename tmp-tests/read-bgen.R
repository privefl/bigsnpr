library(magrittr)
bgen_file <- tempfile(fileext = ".bgen")
system.file("testdata", "bgen_example.rds", package = "bigsnpr") %>%
  readRDS() %>% writeBin(bgen_file, useBytes = TRUE)
system.file("testdata", "bgi_example.rds",  package = "bigsnpr") %>%
  readRDS() %>% writeBin(paste0(bgen_file, ".bgi"), useBytes = TRUE)

header <- readBin(bgen_file, what = 1L, size = 4, n = 5)
as.raw()
rawToChar(as.raw(readBin(bgen_file, what = 1L, size = 1, n = 20)))
rawToChar(as.raw(tail(readBin(bgen_file, what = 1L, size = 1, n = 20), 4)))))

L_H <- header[2]
flags <- rawToBits(tail(readBin(bgen_file, what = raw(), n = L_H + 4), 4))
readBin(flags[1:2], what = integer(), size = 2)
readBin(flags[3:6], what = integer(), size = 4)

N <- header[4]

tmp <- tail(readBin(bgen_file, what = raw(), n = 4 + header[1] + 10 + N), 10 + N)
readBin(tmp, what = integer(), n = 2, size = 4)
unzip()

con <- file(bgen_file, open = "rb")
readBin(con, what = raw(), n = 4 + header[1])
readBin(con, what = raw(), n = readBin(con, what = integer(), size = 2))  # id
readBin(con, what = raw(), n = readBin(con, what = integer(), size = 2))  # rsid
readBin(con, what = raw(), n = readBin(con, what = integer(), size = 2))  # chr

readBin(con, what = integer(), size = 4) # pos
if (readBin(con, what = integer(), size = 2) != 2)
  stop2("Only 2 alleles are allowed.")
readBin(con, what = raw(), n = readBin(con, what = integer(), size = 4))  # a1
readBin(con, what = raw(), n = readBin(con, what = integer(), size = 4))  # a2
readBin(con, what = integer(), size = 4)
D <- readBin(con, what = integer(), size = 4)
close(con)

if (D != (10 + 3 * N))
  stop2("")
