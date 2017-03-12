require(SNPRelate)

root <- "../../Téléchargements/plink_linux_x86_64/celiac300_chr1"
test <- SNPRelate::snpgdsBED2GDS(bed.fn = paste0(root, ".bed"),
                                 fam.fn = paste0(root, ".fam"),
                                 bim.fn = paste0(root, ".bim"),
                                 out.gdsfn = paste0(root, ".gds"))


genofile <- snpgdsOpen(test)

print(system.time(
  GRM <- SNPRelate::snpgdsGRM(genofile)
)) # 15.4 min with 1.7 Gb


require(bigsnpr)

celiac <- snp_attach("backingfiles/celiac300_sub1.rds")

G <- celiac$genotypes
G@code <- c(0:2, NA, 0:2, seq(0, 2, by = 0.01), rep(NA, 48))
print(system.time(
  GRM2 <- big_tcrossprodSelf(G, snp_scaleBinom())
)) # 2.4 min with 5 Gb

all.equal(GRM2$K[1:5, 1:5] / 16325, GRM$grm[1:5, 1:5])

print(system.time(
  hwe <- SNPRelate::snpgdsHWE(genofile)
)) # 2 sec

print(system.time(
  ibs <- SNPRelate::snpgdsIBSNum(genofile)
)) # 2.5 min
ibs$ibs[1:5, 1:5]
