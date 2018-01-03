library(bigsnpr)

gen <- function(n, m) {
  I <- 1:m
  p <- I / (2 * m + 1)

  mat <- outer(I, I, FUN = function(i, j) {
    1 / (abs(i - j) + 1)
  })

  bindata::rmvbin(n, p, bincorr = mat) +
    bindata::rmvbin(n, p, bincorr = mat)
}

set.seed(1)
X <- gen(500, 10)
print(head(X, 20))


plink <- download_plink()
tmp <- tempfile()
fake <- snp_fake(500, 10)
fake$genotypes[] <- X
bed <- snp_writeBed(fake, paste0(tmp, ".bed"))
system(glue::glue("{plink} --bfile {tmp} --out {tmp}",
                  " --indep-pairwise 10 1 0.2"))
