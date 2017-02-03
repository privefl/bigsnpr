require(bigsnpr)

celiac <- snp_attach("backingfiles/celiac300_sub1.bk")
counts <- snp_counts(celiac)
nbNA <- nrow(celiac$genotypes) - colSums(counts$cols.cases + counts$cols.controls)
hist(nbNA, breaks = 100)
table(nbNA)

f <- MASS::fitdistr(nbNA, 'weibull')
nbNA.simu <- round(rweibull(5e5, shape = 0.6, scale = 3))
