test <- snp_readExample()

# generating random phenotypes
test$fam$affection <- sample(c(-1, 1), size = nrow(test$fam),
                         replace = TRUE)

ind <- snp_prune(test)
# ind <- Prune(test, ncores = 2)

R2 <- RsqClass(test$genotypes, test$fam$pheno)

par.save <- par(mfrow = c(1, 2))

plot(R2, type = 'h')
plot(ind, R2[ind], type = "h")

par(par.save)
