require(bigsnpr)
celiac <- snp_attach("backingfiles/celiac300_sub1.rds")
G <- celiac$genotypes
counts <- big_counts(G)
print(counts[4, ])

# all missing values will be 2
G@code <- c(0:2, NA, 0:2, seq(0, 2, by = 0.01), rep(NA, 48))
counts2 <- big_counts(G)
print(counts2[, 1])

print(system.time(
  true <- big_randomSVD(G, big_scale())
)) # 31 sec

print(system.time(
  test <- big_randomSVD(G, big_scale(), ncores = 4)
)) # 20 sec
all.equal(test, true)

print(system.time(
  test3 <- big_randomSVD2(G, big_scale(), ncores = 4)
)) # 21 sec
all.equal(test3, true)



unique(G@code)
rownames(counts2)


X2 <- deepcopy(attach.BM(G), cols = 1:5)
X3 <- X2[,]

table(as.numeric(X3))

nbNA <- big_apply(G, a.FUN = function(x, ind) colSums(is.na(x[, ind])), a.combine = 'c')
all(nbNA == 0)
