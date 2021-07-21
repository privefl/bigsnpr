library(bigsnpr)
NCORES <- nb_cores()

obj.bed <- bed(download_1000G("tmp-data"))
map_1000G <- dplyr::transmute(obj.bed$map,
                              rsid = marker.ID, chr = chromosome, pos = physical.pos,
                              a0 = allele2, a1 = allele1)

dubois <- snp_attach("../Dubois2010_data/FinnuncorrNLITUK1UK3hap300_QC_norel.rds")
map_dubois <- dplyr::transmute(dubois$map,
                               rsid = marker.ID, chr = chromosome, pos = physical.pos,
                               a0 = allele2, a1 = allele1)

matched <- snp_match(cbind(map_dubois, beta = 1), map_1000G, join_by_pos = FALSE)
# ind_match <- vctrs::vec_match(map_dubois[1:2], map_1000G[1:2])
# ind_keep <- which(!is.na(ind_match))

G <- snp_fastImputeSimple(dubois$genotypes, method = "mean2", ncores = NCORES)
set.seed(1); simu <- snp_simuPheno(G, h2 = 0.4, M = 500,
                                   ind.possible = matched[["_NUM_ID_.ss"]])

CHR <- dubois$map$chromosome
POS <- dubois$map$physical.pos
obj.svd <- snp_autoSVD(G, CHR, POS, thr.r2 = 0.1, ncores = NCORES, k = 6)
plot(obj.svd, type = "scores")

gwas <- big_univLinReg(G, simu$pheno, covar.train = obj.svd$u, ncores = NCORES)
plot(gwas)
snp_manhattan(gwas, CHR, POS, npoints = 20e3)

pop <- bigreadr::fread2(paste0(obj.bed$prefix, ".fam2"))[["Super Population"]]
rds <- snp_readBed2(obj.bed$bedfile, ind.row = which(pop == "EUR"),
                    ind.col = matched[["_NUM_ID_"]], backingfile = tempfile())
snp <- snp_attach(rds)
G <- snp$genotypes
CHR <- snp$map$chromosome
POS <- snp$map$physical.pos
snp$map$genetic.dist <- snp_asGeneticPos(CHR, POS, dir = "tmp-data", ncores = NCORES)

ind <- match(simu$set, matched[["_NUM_ID_.ss"]])
causal_beta <- matched$beta[ind] * simu$effects
snp$fam$affection <- scale(G[, ind]) %*% causal_beta

# Verif by compare GWAS
gwas2 <- big_univLinReg(G, snp$fam$affection)
cor(gwas$estim[matched[["_NUM_ID_.ss"]]], gwas2$estim)

sumstats <- cbind.data.frame(
  map_dubois,
  beta = gwas$estim,
  beta_se = gwas$std.err,
  N = length(simu$pheno),
  p = predict(gwas, log10 = FALSE))
sumstats <- sumstats[seq(1, nrow(sumstats), length.out = 50e3), ]
object.size(sumstats) / 1024^2  # 6 MB

bigreadr::fwrite2(sumstats, "tmp-data/public-data3-sumstats.txt")
snp_writeBed(snp, "tmp-data/public-data3.bed",
             ind.col = which(snp$map$marker.ID %in% sumstats$rsid))
zip("data-raw/public-data3.zip",
    paste0("tmp-data/public-data3", c(".bed", ".bim", ".fam", "-sumstats.txt")))
