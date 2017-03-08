require(bigsnpr)

# snp_readBed("../popres/POPRES_data/POPRES_allchr.bed",
#             backingfile = "popres")
popres <- snp_attach("backingfiles/popres.rds")

popresSub <- subset(popres, ind.row = 1:5, ind.col = -1)
popresSub2 <- subset(popres, ind.row = -c(1, 2), ind.col = -1)
