[![Travis-CI Build Status](https://travis-ci.org/privefl/bigsnpr.svg?branch=master)](https://travis-ci.org/privefl/bigsnpr)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/privefl/bigsnpr?branch=master&svg=true)](https://ci.appveyor.com/project/privefl/bigsnpr)
[![Coverage Status](https://img.shields.io/codecov/c/github/privefl/bigsnpr/master.svg)](https://codecov.io/github/privefl/bigsnpr?branch=master)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/bigsnpr)](https://CRAN.R-project.org/package=bigsnpr)
[![DOI](https://zenodo.org/badge/doi/10.1093/bioinformatics/bty185.svg)](http://dx.doi.org/10.1093/bioinformatics/bty185)
 
 
# bigsnpr

{bignspr} is an R package for the analysis of massive SNP arrays. It enhances the features of [package {bigstatsr}](https://privefl.github.io/bigstatsr) for the purpose of analysing genotype data, primarily designed for human genetics.

[Quick demo](https://privefl.github.io/bigsnpr/articles/demo.html)

[**LIST OF FEATURES**](https://privefl.github.io/bigsnpr/reference/index.html)


## Installation

```r
devtools::install_github("privefl/bigsnpr")
```


## Input formats

This package reads *bed*/*bim*/*fam* files (PLINK preferred format) using function `snp_readBed()`. Before reading into this package's special format, quality control and conversion can be done using PLINK, which can be called directly from R using `snp_plinkQC` and `snp_plinkIBDQC`.

This package now also reads **UK Biobank BGEN files** using function `snp_readBGEN()`.

This package uses a class called `bigSNP` for representing SNP data. A `bigSNP` object is just a list with some elements:

- `genotypes`: A [`FBM.code256`](https://privefl.github.io/bigstatsr/reference/FBM.code256-class.html). Rows are samples and columns are SNPs. This stores genotypes calls or dosages (rounded to 2 decimal places).
- `fam`: A `data.frame` containing some information on the SNPs.
- `map`: A `data.frame` giving some information on the individuals.


## Possible upcoming features

- Multiple imputation for GWAS (https://doi.org/10.1371/journal.pgen.1006091).

You can request some feature by opening an issue.


## Bug report

Please open an issue if you find a bug.

If you want help using {bigstatsr}, please post on Stack Overflow with the tag *bigstatsr*.

[How to make a great R reproducible example?](https://stackoverflow.com/q/5963269/6103040)


## Get genes associated with SNPs

```r
devtools::source_gist("42b41d771bbeae63245b8304ef283c70", filename = "get-genes.R")
rsid <- c("rs3934834", "rs3737728", "rs6687776", "rs9651273", "rs4970405",
          "rs12726255", "rs2298217", "rs4970362", "rs9660710", "rs4970420")
snp_gene(rsid)
```
