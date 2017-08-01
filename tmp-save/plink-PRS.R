library(glue)
plink <- "../plink_linux_x86_64/plink"
prefix_bed <- "inst/extdata/example"
write.table2 <- bigsnpr:::write.table2
fread2 <- function(...) data.table::fread(data.table = FALSE, ...)

write.table2(cbind(test$fam[1:2], obj.svd$u), t1 <- tempfile())
system(glue("{plink} --bfile {prefix_bed}",
            " --logistic hide-covar beta",
            " --covar {t1}",
            " --allow-no-sex"))

gwas2 <- fread2("plink.assoc.logistic")
# saveRDS(gwas2$P, "inst/extdata/pval.rds")

system(glue("{plink} --bfile {prefix_bed}",
            " --clump plink.assoc.logistic",
            " --clump-p1 1 --clump-p2 1 --clump-r2 0.2",
            " --allow-no-sex"))

clump <- fread2("plink.clumped")
ind.keep2 <- match(clump$SNP, test$map$marker.ID)
# saveRDS(ind.keep2, "inst/extdata/clumping.rds")

gwas2.keep <- gwas2[ind.keep, ]
write.table2(gwas2.keep, t2 <- tempfile())
write.table2(cbind(paste0("S", seq_along(thrs)), 0, 10^(-thrs)),
             t3 <- tempfile())
system(glue("{plink} --bfile {prefix_bed}",
            " --extract {t2}",
            " --score {t2} 2 4 7 sum",
            " --q-score-range {t3} {t2} 2 9",
            " --allow-no-sex"))

scores <- matrix(NA_real_, nrow(G), length(thrs))
for (i in seq_along(thrs)) {
  scores[, i] <- fread2(glue("plink.S{i}.profile"))$SCORESUM
}
# saveRDS(scores, "inst/extdata/scores-PRS.rds")
