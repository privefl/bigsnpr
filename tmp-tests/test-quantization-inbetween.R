corr0 <- readRDS("../datasets/ldref/LD_with_blocks_chr6.rds")

RSpectra::eigs(corr0, k = 2)$values
# 665.7041 295.2265
RSpectra::eigs(corr0, k = 2, sigma = -1)$values
# -0.1605048 -0.2454353

transform_lookup <- function(x, lookup) {
  ind <- drop(nabor::knn(lookup, x, k = 1)$nn.idx)
  lookup[ind]
}

lookup1 <- seq(-1, 1, length.out = 2^8 - 1)
corr1 <- corr0; corr1@x <- transform_lookup(corr0@x, lookup1)

saveRDS(corr1, "tmp-data/corr_lookup1.rds")
file.size("../datasets/ldref/LD_with_blocks_chr6.rds") /
  file.size("tmp-data/corr_lookup1.rds")  # 9.2

RSpectra::eigs(corr1, k = 2)$values
# 665.6916 295.2457
RSpectra::eigs(corr1, k = 2, sigma = -1)$values
# -0.7826883 -0.8661916


lookup2 <- seq(-1, 1, length.out = 2^16 - 1)
corr2 <- corr0; corr2@x <- transform_lookup(corr0@x, lookup2)

saveRDS(corr2, "tmp-data/corr_lookup2.rds")
file.size("../datasets/ldref/LD_with_blocks_chr6.rds") /
  file.size("tmp-data/corr_lookup2.rds")  # 2.4

RSpectra::eigs(corr2, k = 2)$values
# 665.7041 295.2265
RSpectra::eigs(corr2, k = 2, sigma = -1)$values
# -0.1604970 -0.2454411

FUN <- function(x, power) sign(x) * abs(x)^power
lookup3 <- round(FUN(seq(-1, 1, length.out = 499), 1.25), 4)
corr3 <- corr0; corr3@x <- transform_lookup(corr0@x, lookup3)

saveRDS(corr3, "tmp-data/corr_lookup3.rds")
file.size("../datasets/ldref/LD_with_blocks_chr6.rds") /
  file.size("tmp-data/corr_lookup3.rds")  # 6.8

RSpectra::eigs(corr3, k = 2)$values
# 665.6924 295.2225
RSpectra::eigs(corr3, k = 2, sigma = -1)$values
# -0.2656469 -0.2838187


TRANS <- function(x) sign(x) * floor(127 * abs(x)) / 127
corr4 <- corr0; corr4@x <- TRANS(corr0@x)

saveRDS(corr4, "tmp-data/corr_lookup4.rds")
file.size("../datasets/ldref/LD_with_blocks_chr6.rds") /
  file.size("tmp-data/corr_lookup4.rds")  # 9.2

RSpectra::eigs(corr4, k = 2)$values
# 653.9227 288.8572
RSpectra::eigs(corr4, k = 2, sigma = -1)$values
# -0.9550967 -0.9680745


lookup5 <- round(sinpi(seq(-1, 1, length.out = 399) / 2) * 2**14) / 2**14
cbind(head(lookup5, 10), lookup5[length(lookup5) / 2 + -4:5], tail(lookup5, 10))
corr5 <- corr0; corr5@x <- transform_lookup(corr0@x, lookup5)

saveRDS(corr5, "tmp-data/corr_lookup5.rds")
file.size("../datasets/ldref/LD_with_blocks_chr6.rds") /
  file.size("tmp-data/corr_lookup5.rds")  # 9.7

RSpectra::eigs(corr5, k = 2)$values
# 665.7139 295.2449
RSpectra::eigs(corr5, k = 2, sigma = -1)$values
# -0.7810079 -0.8622470


FUN2 <- Vectorize(function(x) {
  if (x < 0)   return(-Recall(-x))
  if (x > 0.5) return(1 - Recall(1 - x))
  if (x < 0.1) return(-100 * x^3 + 20 * x^2)
  x
})
curve(FUN2(x), n = 1000, from = -1)
lookup6 <- round(FUN2(seq(-1, 1, length.out = 499)) * 2**14) / 2**14
cbind(head(lookup6, 10), lookup6[length(lookup6) / 2 + -4:5], tail(lookup6, 10))
corr6 <- corr0; corr6@x <- transform_lookup(corr0@x, lookup6)

saveRDS(corr6, "tmp-data/corr_lookup6.rds")
file.size("../datasets/ldref/LD_with_blocks_chr6.rds") /
  file.size("tmp-data/corr_lookup6.rds")  # 7.5

RSpectra::eigs(corr6, k = 2)$values
# 665.7054 295.2169
RSpectra::eigs(corr6, k = 2, sigma = -1)$values
# -0.4035851 -0.4270010


lookup7 <- seq(-1, 1, length.out = 999)
corr7 <- corr0; corr7@x <- transform_lookup(corr0@x, lookup7)

saveRDS(corr7, "tmp-data/corr_lookup7.rds")
file.size("../datasets/ldref/LD_with_blocks_chr6.rds") /
  file.size("tmp-data/corr_lookup7.rds")  # 7.3

RSpectra::eigs(corr7, k = 2)$values
# 665.7093 295.2230
RSpectra::eigs(corr7, k = 2, sigma = -1)$values
# -0.3709447 -0.4013577


lookup8 <- round(FUN(seq(-1, 1, length.out = 999), 1.25) * 2**14) / 2**14
cbind(head(lookup8, 10), lookup8[length(lookup8) / 2 + -4:5], tail(lookup8, 10))
corr8 <- corr0; corr8@x <- transform_lookup(corr0@x, lookup8)

saveRDS(corr8, "tmp-data/corr_lookup8.rds")
file.size("../datasets/ldref/LD_with_blocks_chr6.rds") /
  file.size("tmp-data/corr_lookup8.rds")  # 5.9

RSpectra::eigs(corr8, k = 2)$values
# 665.7019 295.2260
RSpectra::eigs(corr8, k = 2, sigma = -1)$values
# -0.1683411 -0.2487258

lookup9 <- round(FUN(seq(-1, 1, length.out = 999), 1.25), 4)
cbind(head(lookup9, 10), lookup9[length(lookup9) / 2 + -4:5], tail(lookup9, 10))
# Precision: 1/2500 close to 0, 1/400 otherwise
corr9 <- corr0; corr9@x <- transform_lookup(corr0@x, lookup9)

saveRDS(corr9, "tmp-data/corr_lookup9.rds")
file.size("../datasets/ldref/LD_with_blocks_chr6.rds") /
  file.size("tmp-data/corr_lookup9.rds")  # 5.8

RSpectra::eigs(corr9, k = 2)$values
# 665.7020 295.2267
RSpectra::eigs(corr9, k = 2, sigma = -1)$values
# -0.1683636 -0.2486049
# PROBABLY THE BEST RIGHT NOW


FUN3 <- Vectorize(function(x) {
  if (x < 0)   return(-Recall(-x))
  (x^1.5 + (1 - (1 - x)^1.5)) / 2
})
curve(FUN3(x), n = 1000, from = -1)
lookup10 <- round(FUN3(seq(-1, 1, length.out = 999)) * 2**14) / 2**14
cbind(head(lookup10, 10), lookup10[length(lookup10) / 2 + -4:5], tail(lookup10, 10))
corr10 <- corr0; corr10@x <- transform_lookup(corr0@x, lookup10)

saveRDS(corr10, "tmp-data/corr_lookup10.rds")
file.size("../datasets/ldref/LD_with_blocks_chr6.rds") /
  file.size("tmp-data/corr_lookup10.rds")  # 6.5

RSpectra::eigs(corr10, k = 2)$values
# 665.7074 295.2285
RSpectra::eigs(corr10, k = 2, sigma = -1)$values
# -0.1732974 -0.2539854


curve(FUN(x, 1.4), n = 1000, from = -1)
lookup11 <- round(FUN(seq(-1, 1, length.out = 999), 1.4) * 2**14) / 2**14
cbind(head(lookup11, 10), lookup11[length(lookup11) / 2 + -4:5], tail(lookup11, 10))
corr11 <- corr0; corr11@x <- transform_lookup(corr0@x, lookup11)

saveRDS(corr11, "tmp-data/corr_lookup11.rds")
file.size("../datasets/ldref/LD_with_blocks_chr6.rds") /
  file.size("tmp-data/corr_lookup11.rds")  # 5.6

RSpectra::eigs(corr11, k = 2)$values
# 665.7037 295.2301
RSpectra::eigs(corr11, k = 2, sigma = -1)$values
# -0.1671062 -0.2501356


lookup12 <- round(FUN(seq(-1, 1, length.out = 749), 1.25), 4)
cbind(head(lookup12, 10), lookup12[length(lookup12) / 2 + -4:5], tail(lookup12, 10))
corr12 <- corr0; corr12@x <- transform_lookup(corr0@x, lookup12)

saveRDS(corr12, "tmp-data/corr_lookup12.rds")
file.size("../datasets/ldref/LD_with_blocks_chr6.rds") /
  file.size("tmp-data/corr_lookup12.rds")  # 6.2

RSpectra::eigs(corr12, k = 2)$values
# 665.7023 295.2255
RSpectra::eigs(corr12, k = 2, sigma = -1)$values
# -0.1755541 -0.2544742


lookup13 <- round(FUN(seq(-1, 1, length.out = 2499), 1.25), 4)
cbind(head(lookup13, 10), lookup13[length(lookup13) / 2 + -4:5], tail(lookup13, 10))
# Precision: at least 1/1000 otherwise, up to 1/10000
corr13 <- corr0; corr13@x <- transform_lookup(corr0@x, lookup13)

saveRDS(corr13, "tmp-data/corr_lookup13.rds")
file.size("../datasets/ldref/LD_with_blocks_chr6.rds") /
  file.size("tmp-data/corr_lookup13.rds")  # 4.7

RSpectra::eigs(corr13, k = 2)$values
# 665.7045 295.2268
RSpectra::eigs(corr13, k = 2, sigma = -1)$values
# -0.1619540 -0.2457191


corr14 <- corr0; corr14@x <- transform_lookup(0.996 * corr0@x, lookup1); diag(corr14) <- 1
sum(corr14@x ==  1) == ncol(corr14)
sum(corr14@x == -1) == 0

saveRDS(corr14, "tmp-data/corr_lookup14.rds")
file.size("../datasets/ldref/LD_with_blocks_chr6.rds") /
  file.size("tmp-data/corr_lookup14.rds")  # 9.2

RSpectra::eigs(corr14, k = 2)$values
# 663.0612 294.0641
RSpectra::eigs(corr14, k = 2, sigma = -1)$values
# -0.8128022 -0.8894983



microbenchmark::microbenchmark(
  readRDS("../datasets/ldref/LD_with_blocks_chr6.rds"),  # 5.2 s
  readRDS("tmp-data/corr_lookup9.rds"),                  # 3.1 s
  readRDS("tmp-data/corr_lookup13.rds"),                 # 3.3 s
  times = 10
)
