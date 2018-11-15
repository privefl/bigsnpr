# http://www.well.ox.ac.uk/~gav/bgen_format/

file <- "tmp-data/first_bytes.bgen"
file_con <- file(file, open = "rb", raw = TRUE)


# The first four bytes
first4 <- readBin(file_con, what = 1L, size = 4)        ## 15,609,036

# The header block
header <- readBin(file_con, what = raw(), n = first4)
writeBin(header, tmp <- tempfile())
header_con <- file(tmp, open = "rb", raw = TRUE)

L_H <- readBin(header_con, what = 1L, size = 4)    ## 20
stopifnot(L_H < first4)

M <- readBin(header_con, what = 1L, size = 4)      ## 7,402,791
N <- readBin(header_con, what = 1L, size = 4)      ## 487,409

magic_number <- readBin(header_con, what = raw(), n = 4)
stopifnot(identical(magic_number, charToRaw("bgen")))

flags <- readBin(header_con, what = raw(), n = 4)
flags_bits <- rawToBits(flags)
flags_bits[1:2]  ## 1 -> use zlib
flags_bits[3:6]  ## 2 -> layout 2
flags_bits[32]   ## 0 -> no sample identifiers
tail(header, 4)


# Variant identifying data
L_id <- readBin(file_con, what = 1L, size = 2, signed = FALSE)    ## 12
id <- readBin(file_con, what = raw(), n = L_id)
rawToChar(id)    ## "1:10177_A_AC"
L_rsid <- readBin(file_con, what = 1L, size = 2, signed = FALSE)  ## 11
rsid <- readBin(file_con, what = raw(), n = L_rsid)
rawToChar(rsid)  ## "rs367896724"
L_chr <- readBin(file_con, what = 1L, size = 2, signed = FALSE)   ## 2
chr <- readBin(file_con, what = raw(), n = L_chr)
rawToChar(chr)   ## "01"
pos <- readBin(file_con, what = 1L, size = 4, signed = FALSE)     ## 10177
K <- readBin(file_con, what = 1L, size = 2, signed = FALSE)       ## 2
stopifnot(K == 2)
L_a1 <- readBin(file_con, what = 1L, size = 4, signed = FALSE)    ## 1
a1 <- readBin(file_con, what = raw(), n = L_a1)
rawToChar(a1)   ## "A"
L_a2 <- readBin(file_con, what = 1L, size = 4, signed = FALSE)    ## 2
a2 <- readBin(file_con, what = raw(), n = L_a2)
rawToChar(a2)   ## "AC"
# possibly more alleles... (but here, only 2!)
C <- readBin(file_con, what = 1L, size = 4, signed = FALSE)       ## 930,596
D <- readBin(file_con, what = 1L, size = 4, signed = FALSE)       ## 1,462,237
probas_compressed <- readBin(file_con, what = raw(), n = C - 4)   ## 78 9c ec ..


# Next variant
L_id <- readBin(file_con, what = 1L, size = 2, signed = FALSE)    ## 12
id <- readBin(file_con, what = raw(), n = L_id)
rawToChar(id)    ## "1:10235_T_TA"


# decompress data
probas <- memDecompress(probas_compressed, type = "gzip")
probas2 <- gzmem::mem_inflate(probas_compressed, format = "zlib", r_guess_size = D + 0)
identical(probas, probas2)
microbenchmark::microbenchmark(
  memDecompress(probas_compressed, type = "gzip"),
  gzmem::mem_inflate(probas_compressed, format = "zlib", r_guess_size = D + 0)
) # absolutely no gain...

# data block
writeBin(probas, tmp <- tempfile())
probas_con <- file(tmp, open = "rb", raw = TRUE)
N2 <- readBin(probas_con, what = 1L, size = 4, signed = FALSE)
stopifnot(N2 == N)
K2 <- readBin(probas_con, what = 1L, size = 2, signed = FALSE)
stopifnot(K2 == K)
P_min <- readBin(probas_con, what = 1L, size = 1, signed = FALSE)  ## 2
P_max <- readBin(probas_con, what = 1L, size = 1, signed = FALSE)  ## 2
ploidy <- readBin(probas_con, what = 1L, size = 1, signed = FALSE, n = N)
stopifnot(all(ploidy == 2))  ## diploid + no missing value
phased <- readBin(probas_con, what = 1L, size = 1, signed = FALSE) ## no
B <- readBin(probas_con, what = 1L, size = 1, signed = FALSE)      ## 8
stopifnot(B == 8)
probas_read <- readBin(probas_con, what = 1L, size = 1, signed = FALSE, n = 5 * N)
stopifnot(length(probas_read) == 2 * N)
probas_read
summary(dosage <- 2 - (2 * probas_read[c(TRUE, FALSE)] + probas_read[c(FALSE, TRUE)]) / 255)
stopifnot(D == (4 + 2 + 1 + 1 + N + 1 + 1 + 2 * N))
stopifnot(D == (10 + 3 * N))
