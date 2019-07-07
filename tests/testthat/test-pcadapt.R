################################################################################

context("PCADAPT") # Need to change tolerances

################################################################################

bigsnp <- snp_attachExtdata()
expect_s3_class(bigsnp, "bigSNP")
G <- bigsnp$genotypes
expect_s4_class(G, "FBM.code256")

################################################################################

test_that("Same as pcadapt", {

  # skip_on_cran()

  ################################################################################

  tmpfile <- tempfile(fileext = ".pcadapt")
  write.table(t(G[]), tmpfile, quote = FALSE, sep = " ",
              row.names = FALSE, col.names = FALSE)

  bed <- pcadapt::read.pcadapt(tmpfile, type = "pcadapt")
  obj.pcadapt <- pcadapt::pcadapt(bed, K = 10, min.maf = 0)

  ################################################################################

  obj.svd <- big_SVD(G, snp_scaleBinom())
  test <- bigsnpr:::multLinReg(G, rows_along(G), cols_along(G), obj.svd$u)

  expect_equal(obj.svd$center / 2, obj.pcadapt$af)
  expect_equal(obj.svd$d, obj.pcadapt$singular.values * sqrt((nrow(G) - 1) * ncol(G)))
  expect_equal(obj.svd$u, obj.pcadapt$scores,   tolerance = 1e-4)
  expect_equal(obj.svd$v, obj.pcadapt$loadings, tolerance = 1e-2)
  expect_equal(test,      obj.pcadapt$zscores,  tolerance = 1e-1)
  # expect_equal(abs(cov(obj.svd$u, obj.pcadapt$scores)), diag(10) / (nrow(G) - 1))
  # expect_equal(abs(cor(obj.svd$v, obj.pcadapt$loadings)), diag(10),
  #              tolerance = 1e-2, check.attributes = FALSE)
  # expect_equal(abs(cor(test, obj.pcadapt$zscores)), diag(10),
  #              tolerance = 1e-2, check.attributes = FALSE)

  ################################################################################

  obj.gwas <- snp_pcadapt(G, obj.svd$u, ncores = sample(1:2, 1))
  expect_equal(bigsnpr:::getLambdaGC(obj.gwas), obj.pcadapt$gif,
               tolerance = 1e-2)
  expect_equal(snp_gc(obj.gwas)[[1]], as.numeric(obj.pcadapt$stat),
               tolerance = 1e-1)
  expect_equal(predict(snp_gc(obj.gwas), log10 = FALSE), obj.pcadapt$pvalues,
               tolerance = 1e-1)

  ################################################################################

  expect_s3_class(snp_qq(obj.gwas, lambdaGC = FALSE), "ggplot")
  p <- snp_qq(snp_gc(obj.gwas))
  expect_equal(as.character(p$labels$subtitle), "lambda[GC] == 1")
  expect_s3_class(p, "ggplot")
  p2 <- snp_manhattan(obj.gwas, infos.chr = bigsnp$map$chromosome,
                      infos.pos = bigsnp$map$physical.pos,
                      npoints = 2000)
  expect_s3_class(p2, "ggplot")

  ################################################################################

  # K = 1
  obj.pcadapt <- pcadapt::pcadapt(bed, K = 1, min.maf = 0)
  obj.gwas    <- snp_pcadapt(G, obj.svd$u[, 1], ncores = sample(1:2, 1))
  obj.gwas.gc <- snp_gc(obj.gwas)

  snp_qq(obj.gwas)
  snp_qq(obj.gwas.gc)

  tmp <- obj.pcadapt$scores[, 1]; names(tmp) <- NULL
  expect_equal(tmp, obj.svd$u[, 1], tolerance = 1e-6)
  expect_equal(as.vector(obj.pcadapt$zscores), obj.gwas$score, tolerance = 1e-2)
  # plot(obj.pcadapt$pvalues, predict(obj.gwas.gc, log10 = FALSE))

  # expect_equal(bigsnpr:::getLambdaGC(obj.gwas), obj.pcadapt$gif,
  #              tolerance = 1e-2)
  # expect_equal(snp_gc(obj.gwas)[[1]], as.numeric(obj.pcadapt$stat),
  #              tolerance = 1e-2)
  # expect_equal(predict(snp_gc(obj.gwas), log10 = FALSE), obj.pcadapt$pvalues,
  #              tolerance = 1e-2)

  ################################################################################

})

