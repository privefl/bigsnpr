## bigsnpr 0.10.2

- Parameter `is.size.in.bp` is deprecated.

## bigsnpr 0.10.1

- Add parameter `read_as` for `snp_readBGEN()`. It is now possible to sample BGEN probabilities as random hard calls using `read_as = "random"`. Default remains reading probabilities as dosages.

## bigsnpr 0.10.0

- For memory-mapping, now use *mio* instead of *boost*.

- `snp_clumping()` (and `snp_autoSVD()`) now has a `size` that is inversely proportional to `thr.r2`.

- `snp_pruning()` is deprecated (and will be removed someday); now always use `snp_clumping()`.

## bigsnpr 0.9.0

- When reading bed files, switch reading of Os and 2s to be consistent with other software.

## bigsnpr 0.8.2

- Add function `snp_assocBGEN()` for computing quick association tests from BGEN files. Could be useful for quick screening of useful SNPs to read in bigSNP format. This function might be improved in the future.

## bigsnpr 0.8.1

- Change url to download PLINK 1.9.

## bigsnpr 0.8.0

- Add function `snp_readBGEN()` to read UK Biobank BGEN files in `bigSNP` format.

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

- Faster defaults + possibility to estimate correlations based on a subset of individuals for `snp_fastImpute`. Also store information in an FBM (instead of a data frame) so that imputation can be done by parts (you can stop the imputation by killing the R processes and come back to it later). Note that the defaults used in the *Bioinformatics* paper were `alpha = 0.02` and `size = 500` (instead of `1e-4` and `200` now, respectively). These new defaults are more stringent on the SNPs that are used, which makes the imputation faster (30 min instead of 42-48 min), without impacting accuracy (still 4.7-4.8% of errors).

## bigsnpr 0.2.5

- **This package won't be on CRAN**.

## bigsnpr 0.2.4

- No longer download PLINK automatically (because it is a CRAN policy violation).
