require(SNPRelate)

# root <- "../../Téléchargements/plink_linux_x86_64/celiac300_chr1"
# root <- "../../plink_linux_x86_64/celiac300_sub"
# test <- SNPRelate::snpgdsBED2GDS(bed.fn = paste0(root, ".bed"),
#                                  fam.fn = paste0(root, ".fam"),
#                                  bim.fn = paste0(root, ".bim"),
#                                  out.gdsfn = paste0(root, ".gds"))


genofile <- snpgdsOpen(paste0(root, ".gds"))

print(system.time(
  relate <- snpgdsIBDMoM(genofile, num.thread = 6)
)) # 7.5 min (6 cores) -> 12 min with PLINK (11 cores)

genome <- data.table::fread("../../plink_linux_x86_64/plink.genome")
require(dplyr)
ind <- as.matrix(transmute(genome, IID1 = match(IID1, relate$sample.id),
                           IID2 = match(IID2, relate$sample.id)))

all.equal(relate$k0[ind], genome$Z0) # Mean relative difference: 3.424118e-05
all.equal(relate$k1[ind], genome$Z1) # Mean relative difference: 0.0001229569

kinship <- (1 - relate$k0[ind] - relate$k1[ind] / 2) / 2
all.equal(2 * kinship, genome$PI_HAT) # Mean relative difference: 0.0002069291
all.equal(rowSums(genome[, paste0("Z", 0:2)]), rep(1, 237)) # Mean relative difference: 8.444332e-05

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
ibs$ibs0[1:5, 1:5]
