################################################################################

context("PCADAPT") # Need to change tolerances

################################################################################

bigsnp <- snp_attachExtdata()
expect_s3_class(bigsnp, "bigSNP")
G <- bigsnp$genotypes
expect_s4_class(G, "big.matrix.descriptor")
expect_s4_class(G, "BM.code.descriptor")

################################################################################

tmpfile <- tempfile(fileext = ".pcadapt")
write.table(t(attach.BM(G)[]), tmpfile, quote = FALSE, sep = " ",
            row.names = FALSE, col.names = FALSE)

obj.pcadapt <- pcadapt::pcadapt(tmpfile, K = 10, min.maf = 0)

################################################################################

obj.svd <- big_SVD(G, snp_scaleBinom())
test <- bigsnpr:::linRegPcadapt_cpp(attach.BM(G), obj.svd$u,
                                    rows_along(G), cols_along(G))

expect_equal(pmin(obj.svd$means, 2 - obj.svd$means) / 2, obj.pcadapt$maf,
             tolerance = 1e-6)
expect_equal(abs(cor(obj.svd$v, obj.pcadapt$loadings)), diag(10),
             tolerance = 1e-1, check.attributes = FALSE)
expect_equal(abs(cov(obj.svd$u, obj.pcadapt$scores)), diag(10) / (nrow(G) - 1))
expect_equal(abs(cor(test, obj.pcadapt$zscores)), diag(10),
             tolerance = 1e-1, check.attributes = FALSE)

################################################################################

obj.gwas <- snp_pcadapt(G, obj.svd$u)
expect_equal(bigsnpr:::getLambdaGC(obj.gwas), obj.pcadapt$gif,
             tolerance = 1e-2)
expect_equal(snp_gc(obj.gwas)[[1]], as.numeric(obj.pcadapt$stat),
             tolerance = 1e-2)
expect_equal(predict(snp_gc(obj.gwas), log10 = FALSE), obj.pcadapt$pvalues,
             tolerance = 1e-2)

################################################################################
