PKGs <- c("bigmemory", "biglasso", "bigstatsr", "bigsnpr")
pacman::p_load_current_gh(paste("privefl", PKGs, sep = "/"))


celiac <- snp_attach("backingfiles/celiac.bk")
counts <- snp_counts(celiac, FALSE)
nbNA <- nrow(celiac$genotypes) - colSums(counts$cols)
hist(nbNA, breaks = 100)
table(nbNA)

nbNA2 <- nbNA[nbNA < 1000]

f <- MASS::fitdistr(nbNA2+0.01, 'weibull')
nbNA.simu <- round(rweibull(5e5, shape = 0.6, scale = 3))

require(VGAM)
require(fitdistrplus)
f2 <- fitdist(nbNA2, "betabinom.ab", start = list(shape1 = 1, shape2 = 1),
              fix.arg = list(size = nrow(celiac$genotypes)), discrete = TRUE)

nbNA.simu2 <- rbetabinom.ab(5e5, size = nrow(celiac$genotypes),
                           shape1 = f2$estimate[1], shape2 = f2$estimate[2])
table(nbNA.simu2)
