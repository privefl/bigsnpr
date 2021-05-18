<!-- badges: start -->
[![R build status](https://github.com/privefl/bigsnpr/workflows/R-CMD-check/badge.svg)](https://github.com/privefl/bigsnpr/actions)
[![Travis-CI Build Status](https://travis-ci.org/privefl/bigsnpr.svg?branch=master)](https://travis-ci.org/privefl/bigsnpr)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/privefl/bigsnpr?branch=master&svg=true)](https://ci.appveyor.com/project/privefl/bigsnpr)
[![Coverage Status](https://img.shields.io/codecov/c/github/privefl/bigsnpr/master.svg)](https://codecov.io/github/privefl/bigsnpr)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/bigsnpr)](https://CRAN.R-project.org/package=bigsnpr)
[![DOI](https://zenodo.org/badge/doi/10.1093/bioinformatics/bty185.svg)](http://dx.doi.org/10.1093/bioinformatics/bty185)
<!-- badges: end -->
 
 
# bigsnpr

{bigsnpr} is an R package for the analysis of massive SNP arrays, primarily designed for human genetics. It enhances the features of [package {bigstatsr}](https://privefl.github.io/bigstatsr/) for the purpose of analyzing genotype data.

[Quick demo](https://privefl.github.io/bigsnpr/articles/demo.html)

[**LIST OF FEATURES**](https://privefl.github.io/bigsnpr/reference/index.html)

**Note that most of the algorithms of this package don't handle missing values.** You can use `snp_fastImpute()` (taking a few hours for a chip of 15K x 300K) and `snp_fastImputeSimple()` (taking a few minutes only) to impute missing values of *genotyped* variants.

**New!** Package {bigsnpr} now provides functions that directly work on bed files with a few missing values. See new paper "Efficient toolkit implementing..".


## Installation

In R, run

```r
# install.packages("remotes")
remotes::install_github("privefl/bigsnpr")
```

or for the CRAN version

```r
install.packages("bigsnpr")
```


## Input formats

This package reads *bed*/*bim*/*fam* files (PLINK preferred format) using functions `snp_readBed()` and `snp_readBed2()`. Before reading into this package's special format, quality control and conversion can be done using PLINK, which can be called directly from R using `snp_plinkQC()` and `snp_plinkKINGQC()`.

This package can also read **UK Biobank BGEN files** using function `snp_readBGEN()`. This function takes around 40 minutes to read 1M variants for 400K individuals using 15 cores.

This package uses a class called `bigSNP` for representing SNP data. A `bigSNP` object is a list with some elements:

- `genotypes`: A [`FBM.code256`](https://privefl.github.io/bigstatsr/reference/FBM.code256-class.html). Rows are samples and columns are variants. This stores genotype calls or **dosages** (rounded to 2 decimal places).
- `fam`: A `data.frame` with some information on the individuals.
- `map`: A `data.frame` with some information on the variants.

**New!** Package {bigsnpr} now provides functions that directly work on bed files with a few missing values. See new paper "Efficient toolkit implementing..".


## Polygenic scores

Polygenic scores are one of the main focus of this package. There are 3 main methods currently available:

- Penalized regressions with individual-level data (see [paper](https://doi.org/10.1534/genetics.119.302019) and [tutorial](https://privefl.github.io/bigstatsr/articles/penalized-regressions.html))

- Clumping and Thresholding (C+T) and Stacked C+T (SCT) with summary statistics and individual level data (see [paper](https://doi.org/10.1016/j.ajhg.2019.11.001) and [tutorial](https://privefl.github.io/bigsnpr/articles/SCT.html)).

- LDpred2 with summary statistics (see [preprint](https://doi.org/10.1101/2020.04.28.066720) and [tutorial](https://privefl.github.io/bigsnpr/articles/LDpred2.html))


## Possible upcoming features

- Multiple imputation for GWAS (https://doi.org/10.1371/journal.pgen.1006091).

- More interactive (visual) QC.

You can request some feature by opening an issue.


## Bug report

[How to make a great R reproducible example?](https://stackoverflow.com/q/5963269/6103040)

Please open an issue if you find a bug.

If you want help using {bigstatsr}, please open an issue on [{bigstatsr}'s repo](https://github.com/privefl/bigstatsr/issues) or post on Stack Overflow with the tag *bigstatsr*.

I will always redirect you to GitHub issues if you email me, so that others can benefit from our discussion.


## References

- Privé, Florian, et al. ["Efficient analysis of large-scale genome-wide data with two R packages: bigstatsr and bigsnpr."](https://doi.org/10.1093/bioinformatics/bty185) Bioinformatics 34.16 (2018): 2781-2787.

- Privé, Florian, et al. ["Efficient implementation of penalized regression for genetic risk prediction."](https://doi.org/10.1534/genetics.119.302019) Genetics 212.1 (2019): 65-74.

- Privé, Florian, et al. ["Making the most of Clumping and Thresholding for polygenic scores."](https://doi.org/10.1016/j.ajhg.2019.11.001) The American Journal of Human Genetics 105.6 (2019): 1213-1221.

- Privé, Florian, et al. ["Efficient toolkit implementing best practices for principal component analysis of population genetic data."](https://doi.org/10.1093/bioinformatics/btaa520) Bioinformatics (2020).

- Privé, Florian, et al. ["LDpred2: better, faster, stronger."](https://doi.org/10.1101/2020.04.28.066720) BioRxiv (2020).

<br>
