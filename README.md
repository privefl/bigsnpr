[![Travis-CI Build Status](https://travis-ci.org/privefl/bigsnpr.svg?branch=master)](https://travis-ci.org/privefl/bigsnpr)
[![Coverage Status](https://img.shields.io/codecov/c/github/privefl/bigsnpr/master.svg)](https://codecov.io/github/privefl/bigsnpr?branch=master)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/bigsnpr)](http://cran.r-project.org/package=bigsnpr)
 
 
# bigsnpr

**bignspr** is an R Package for the analysis of massive SNP arrays. It enhances the features of [package **bigstatsr**](https://privefl.github.io/bigstatsr) for the purpose of analysing genotype data.

[**LIST OF FEATURES**](https://privefl.github.io/bigsnpr/reference/index.html)


## This package is in beta testing

Any bug report is welcomed.


## Installation

For now, you can install this package using

```r
devtools::install_github("privefl/bigsnpr")
```


## Input format

For now, this package only read *bed*/*bim*/*fam* files (PLINK preferred format) using `snp_readBed`. Before reading into this package's special format, quality control and conversion can be done using PLINK, which can be called directly from R using `snp_plinkQC` and `snp_plinkIBDQC`.

I use a class called `bigSNP` for representing infos on massive SNP arrays. One `bigSNP` has at least 3 elements:
- `genotypes`: A `BM.code.descriptor` (see [package **bigstastr**](https://privefl.github.io/bigsnpr/reference/bigSNP-class.html)) which describes a special `big.matrix` encoded with type `raw` (one byte unsigned integer), representing genotype calls and possibly imputed allele dosages. Rows are samples and columns are SNPs.
- `fam`: A `data.frame` containing some information on the SNPs (read from the ".fam" file).
- `map`: A `data.frame` giving some information on the individuals (read from a ".bim" file).


## Possible upcoming features

- Support for other input formats. Note that there is room for coding **allele dosages** (rounded to two decimal places). See `bigsnpr:::CODE_DOSAGE`.
- Fast imputation algorithm which doesn't require reference panels.
- Imputation probabilities and multiple imputation.
- An interactive QC procedure (call rates, difference of missingness between cases and controls, MAF cutoff, relatedness, HWE, autosomal only, others?). 
- proper integration of haploid species.


## Code of conduct

Please note that this project is released with a [Contributor Code of Conduct](https://github.com/privefl/bigsnpr/blob/master/code_of_conduct.md). 
By participating in this project you agree to abide by its terms.
