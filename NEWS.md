## bigsnpr 0.7.0

- Add parameter `is.size.in.bp` to `snp_autoSVD()` for the clumping part.

- Change the threshold of outlier detection in `snp_autoSVD()` (it now detects less outliers). See the documentation details if you don't have any information about SNPs.

## bigsnpr 0.6

- Keep up with {bigstatsr}.

## bigsnpr 0.3.1

- Provide function `snp_gene` (as a gist) to get genes corresponding to 'rs' SNP IDs thanks to package {rsnps} from rOpenSci. See README.

## bigsnpr 0.3.0

- **Package {bigsnpr} is published in [Bioinformatics](http://dx.doi.org/10.1093/bioinformatics/bty185)**.

## bigsnpr 0.2.7

- Faster defaults + possibility to estimate correlations based on a subset of individuals for `snp_fastImpute`. Also store infos in an FBM (instead of a data frame) so that imputation can be done by parts (you can stop the imputation by killing the R processes and come back to it later). Note that the defaults used in the *Bioinformatics* paper were `alpha = 0.02` and `size = 500` (instead of `1e-4` and `200` now, respectively). These new defaults are more stringent on the SNPs that are used, which makes the imputation faster (30 min instead of 42-48 min), without impacting accuracy (still 4.7-4.8% of errors).

## bigsnpr 0.2.5

- **This package won't be on CRAN**.

## bigsnpr 0.2.4

- No longer download PLINK automatically (because it is a CRAN policy violation).
