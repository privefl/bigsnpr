celiac300 <- snp_attach("backingfiles/celiac300.rds")
chrs <- rle(celiac300$map$chromosome)
chrs$values
chrs$lengths

hh <- data.table::fread("../../Téléchargements/plink_linux_x86_64/plink.hh",
                       data.table = FALSE)
ind <- match(unique(hh$V3), celiac300$map$marker.ID)
celiac300$map$chromosome[ind]
X <- attach.BM(celiac300$genotypes)
X[, ind]
counts.ca <- big_counts(X, ind.row = which(celiac300$fam$affection == 2),
                        ind.col = ind)
counts.co <- big_counts(X, ind.row = which(celiac300$fam$affection == 1),
                        ind.col = ind)


counts.ca[, 1:5]
counts.co[, 1:5]

ca.fe <- (celiac300$fam$affection == 2) & (celiac300$fam$sex == 2)
counts.ca.fe <- big_counts(X, ind.row = which(ca.fe), ind.col = ind)
counts.ca.fe[, 1:10]

ca.ma <- (celiac300$fam$affection == 2) & (celiac300$fam$sex == 1)
counts.ca.ma <- big_counts(X, ind.row = which(ca.ma), ind.col = ind)
counts.ca.ma[, 1:10]

co.fe <- (celiac300$fam$affection == 1) & (celiac300$fam$sex == 2)
counts.co.fe <- big_counts(X, ind.row = which(co.fe), ind.col = ind)
counts.co.fe[, 1:10]

co.ma <- (celiac300$fam$affection == 1) & (celiac300$fam$sex == 1)
counts.co.ma <- big_counts(X, ind.row = which(co.ma), ind.col = ind)
counts.co.ma[, 1:10]

counts.ma <- big_counts(X, ind.row = which(celiac300$fam$sex == 1), ind.col = ind)
sum(counts.ma["1", ])
