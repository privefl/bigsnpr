library(bigsnpr)
celiac <- snp_attach("../paper2-PRS/backingfiles/celiacQC.rds")
G <- celiac$genotypes
y <- celiac$fam$affection - 1L
CHR <- celiac$map$chromosome
POS <- celiac$map$physical.pos

ind.row <- rows_along(G)
ind.chr <- which(CHR == 1)
S <- NULL
infos.pos <- POS
size <- 5000
thr.r2 <- 0.95

# cache some computations
stats <- big_colstats(G, ind.row = ind.row, ind.col = ind.chr)
n <- length(ind.row)
denoX <- (n - 1) * stats$var

# statistic to prioritize SNPs
if (is.null(S)) {
  af <- stats$sum / (2 * n)
  S.chr <- pmin(af, 1 - af)
} else {
  S.chr <- S[ind.chr]
}
ord.chr <- order(S.chr, decreasing = TRUE)

# main algo
system.time(
  keep1 <- clumping_chr(
    G,
    rowInd = ind.row,
    colInd = ind.chr,
    ordInd = ord.chr,
    pos    = `if`(is.null(infos.pos), 1000L * seq_along(ind.chr),
                  infos.pos[ind.chr]),
    sumX   = stats$sum,
    denoX  = denoX,
    size   = size * 1000L, # in bp
    thr    = thr.r2
  )
) # 231


library(Matrix)
M <- length(ind.chr)
spcor <- sparseMatrix(i = integer(), j = integer(), x = double(), dims = c(M, M))
system.time(
  res <- clumping_chr_cached(
    G, spcor,
    rowInd = ind.row,
    colInd = ind.chr,
    ordInd = ord.chr,
    pos    = `if`(is.null(infos.pos), 1000L * seq_along(ind.chr),
                  infos.pos[ind.chr]),
    sumX   = stats$sum,
    denoX  = denoX,
    size   = size * 1000L, # in bp
    thr    = thr.r2
  )
) # 239 -> 6
identical(res[[1]], keep1)
spcor <- res[[2]]
spcor[1:10, 1:6]
