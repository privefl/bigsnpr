[![Travis-CI Build Status](https://travis-ci.org/privefl/bigsnpr.svg?branch=master)](https://travis-ci.org/privefl/bigsnpr)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/privefl/bigsnpr?branch=master&svg=true)](https://ci.appveyor.com/project/privefl/bigsnpr)
[![Coverage Status](https://img.shields.io/codecov/c/github/privefl/bigsnpr/master.svg)](https://codecov.io/github/privefl/bigsnpr?branch=master)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/bigsnpr)](https://CRAN.R-project.org/package=bigsnpr)
 
 
# bigsnpr

**bignspr** is an R Package for the analysis of massive SNP arrays. It enhances the features of [package **bigstatsr**](https://privefl.github.io/bigstatsr) for the purpose of analysing genotype data.

[Quick demo](https://privefl.github.io/bigsnpr/articles/demo.html)

[**LIST OF FEATURES**](https://privefl.github.io/bigsnpr/reference/index.html)


## Installation

```r
# Not on CRAN for now because of download_plink()
# For the current version
devtools::install_github("privefl/bigsnpr")
```


## Input format

For now, this package only read *bed*/*bim*/*fam* files (PLINK preferred format) using `snp_readBed`. Before reading into this package's special format, quality control and conversion can be done using PLINK, which can be called directly from R using `snp_plinkQC` and `snp_plinkIBDQC`.

I use a class called `bigSNP` for representing infos on massive SNP arrays. One `bigSNP` has at least 3 elements:
- `genotypes`: A [`FBM.code256`](https://privefl.github.io/bigstatsr/reference/FBM.code256-class.html). Rows are samples and columns are SNPs. This corresponds to the "bed" file, but each element is encoded on 8 bits rather than only 2 bits for PLINK binary files, which allows for storing more information, without taking too much disk space.
- `fam`: A `data.frame` containing some information on the SNPs (read from the ".fam" file).
- `map`: A `data.frame` giving some information on the individuals (read from the ".bim" file).


## Possible upcoming features

- Support for other input formats. Note that there is room for coding **allele dosages** (rounded to two decimal places). See [this vignette](https://privefl.github.io/bigsnpr/articles/dosage.html).
- Imputation of probabilities and multiple imputation.
- An interactive QC procedure (call rates, difference of missingness between cases and controls, MAF cutoff, relatedness, HWE, autosomal only, others?). 
- Proper integration of haploid species.


## Bug report

Please open an issue if you find a bug.
If you want help using **bigstatsr**, please post on Stack Overflow with the tag *bigstatsr* (not yet created). [How to make a great R reproducible example?](https://stackoverflow.com/q/5963269/6103040)


## Code of conduct

Please note that this project is released with a [Contributor Code of Conduct](https://github.com/privefl/bigsnpr/blob/master/code_of_conduct.md). 
By participating in this project you agree to abide by its terms.
