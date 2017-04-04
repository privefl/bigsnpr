[![Travis-CI Build Status](https://travis-ci.org/privefl/bigsnpr.svg?branch=master)](https://travis-ci.org/privefl/bigsnpr)
[![Coverage Status](https://img.shields.io/codecov/c/github/privefl/bigsnpr/master.svg)](https://codecov.io/github/privefl/bigsnpr?branch=master)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/bigsnpr)](http://cran.r-project.org/package=bigsnpr)
 
# bigsnpr

bignspr is an R Package for the analysis of massive SNP arrays.

## This package is under heavy dev

## Input format

For now, this package only read *bed*/*bim*/*fam* files (PLINK preferred format) using `snp_readBed`. Before reading into this package's special format, quality control and conversion can be done using PLINK, which can be called directly from R using `snp_plinkQC` and `snp_plinkIBDQC`.

#TODO: describe bigSNP

## Possible upcoming features

- Support for other input formats. Note that there is room for coding allele dosages (rounded to two decimal places). See `bigsnpr:::CODE_DOSAGE`.
- Fast imputation algorithm which doesn't require reference panels.
- Imputation probabilities and multiple imputation.
- An interactive QC procedure (call rates, difference of missingness between cases and controls, MAF cutoff, relatedness, HWE, autosomal only, others?). 
- proper integration of haploid species.
