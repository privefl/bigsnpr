## bigsnpr 1.12.14

- Add a new `min.maf = 0.02` parameter to `snp_autoSVD()` and `bed_autoSVD()`. Then variants are now discarded when they have either a small MAC or a small MAF.

## bigsnpr 1.12.13

- Can now use matrix accessors for class `bed_light` as well.

## bigsnpr 1.12.12

- Minor improvements to `snp_autoSVD()` and `bed_autoSVD()`: 
    - error when `min.mac = 0`,
    - return a better `attr(, "lrldr")`.

## bigsnpr 1.12.11

- In functions `snp_autoSVD()` and `bed_autoSVD()`, now perform the MAC thresholding before the clumping step. This reordering should not change results, but this should be faster now.

## bigsnpr 1.12.10

- In function `snp_ancestry_summary()`, add parameter `sum_to_one` to optionally allows for ancestry coefficients to have a sum lower than 1 (when `FALSE`; default is `TRUE`).

## bigsnpr 1.12.9

- In function `snp_modifyBuild()`, you can now provide `local_chain` as a vector of two, for when using `check_reverse`. You can now also modify the `base_url` from where to download the chain files.

## bigsnpr 1.12.8

- In function `snp_ancestry_summary()`, now also report correlations between input frequencies and each reference frequencies as well as predicted frequencies. Also add a new parameter `min_cor` to error when the latter correlation is too small.

## bigsnpr 1.12.7

- In function `snp_modifyBuild()`, fix a ftp broken link, and add the possibility to use a local chain file specified by the new parameter `local_chain`.

## bigsnpr 1.12.6

- Fix issue with `snp_subset()` when either `$fam` or `$map` are missing.

## bigsnpr 1.12.5

- Add function `snp_asGeneticPos2()` (and `download_genetic_map()`) where you can provide any reference genetic map as a data frame. This function uses linear interpolation to transform physical positions (in bp) to genetic positions (in cM).

## bigsnpr 1.12.4

- Add parameter `p_bounds` in LDpred2-auto to provide bounds for the estimate of the polygenicity p.

## bigsnpr 1.12.3

- Fix sampling issue of `snp_simuPheno()` when `length(ind.possible)` is 1.

## bigsnpr 1.12.2

- Implement matrix accessors `[,]` for bed objects.

## bigsnpr 1.12.1

- Use a safer detection of strong divergences in LDpred2 and lassosum2.

## bigsnpr 1.12.0

- Add new parameter `ind.corr` to `snp_lassosum2()`, `snp_ldpred2_grid()` and `snp_ldpred2_auto()` to be able to use a subset of `corr` without making a copy of it.

- Add new parameter `ind.beta` to `snp_ldsc2()` to use a subset of the full LD scores corresponding to `df_beta`.

## bigsnpr 1.11.12

- Add new parameter `pos_scaled` to `snp_ldsplit()`.

## bigsnpr 1.11.11

- Fix C++ code that used integers to store the positions for clumping.

## bigsnpr 1.11.10

- Add other architectures (AMD / ARM) as options for PLINK2.

## bigsnpr 1.11.9

- Add option `use_MLE` in LDpred2-auto to allow, when using `FALSE`, for running LDpred2-auto as in previous versions (e.g. v1.10.8), which did not include alpha in the model. Default is `TRUE`.

## bigsnpr 1.11.8

- Detect strong divergence in LDpred2-auto, and return missing values in that case.

## bigsnpr 1.11.7

- Fix font rendering issue of `>=` in subtitle of `snp_manhattan()`.

## bigsnpr 1.11.5

- Autocomplete PLINK builds to be downloaded (fix #383).

## bigsnpr 1.11.4

- Extend and improve LDpred2-auto to allow for estimating $h^2$, $p$, and $\alpha$, a new third parameter that controls how expected effect sizes relate to minor allele frequencies.

## bigsnpr 1.11.3

- Can now run `snp_ldsc2()` with `corr` as an SFBM.

## bigsnpr 1.11.2

- Now use a sparse format for sampling betas returned in LDpred2-auto, instead of a dense matrix that could require quite some memory to store.

## bigsnpr 1.11.0

- Add two new parameters to `snp_simuPheno()`: `alpha` and `prob`.

## bigsnpr 1.10.7

- Fix a liftOver error in `snp_modifyBuild()`.

## bigsnpr 1.10.6

- Better `snp_ldsplit()`:
    - also return `$cost2`, the sum of squared sizes of the blocks,
    - for equivalent splits (with the same cost), now return the one that also minimizes cost2,
    - now return unique splits only (e.g. could get equivalent splits with different `max_size`).

## bigsnpr 1.10.5

- Slightly change the default parameters of lassosum2: 
    - `delta` from `c(0.001, 0.005, 0.02, 0.1, 0.6, 3)` to `c(0.001, 0.01, 0.1, 1)`,
    - `nlambda` from 20 to 30,
    - `maxiter` from 500 to 1000.
    
- Add a penalty multiplicative factor for delta and lambda to regularize variants with smaller GWAS sample sizes more (when they are different, as in meta-analyses with different sets of variants).

## bigsnpr 1.10.4

- Now use the same updating strategy for residuals in LDpred2 as in lassosum2. This can make LDpred2-grid and LDpred2-auto an order of magnitude faster, especially for small p.

## bigsnpr 1.10.2

- Better `snp_modifyBuild()`: more variants should be mapped + add some QC on the mapping (a position is not mapped to more than one, the chromosome is the same, and possibly check whether we can go back to the initial position -> cf. https://doi.org/10.1093/nargab/lqaa054).

## bigsnpr 1.10.1

- Add two new parameters to `snp_ldsplit()`: `max_r2`, the maximum squared correlation allowed outside blocks, and `max_cost`, the maximum cost of reported solutions (i.e. the sum of all squared correlations outside blocks). Using `max_r2` offers an extra guarantee that the splitting is very good, and makes the function much faster by discarding lots of possible splits. 

## bigsnpr 1.10.0

- LDpred2-grid does not use OpenMP for parallelism anymore, it now simply uses multiple R processes. 

- LDpred2-grid and LDpred2-auto can now make use of `set.seed()` to get reproducible results. Note that LDpred2-inf and lassosum2 do not use any sampling.

## bigsnpr 1.9.13

- Enforce `scipen = 50` when writing files to turn off scientific format (e.g. for physical positions stored as `double`).

## bigsnpr 1.9.12 & bigsparser 0.6

- Use a better strategy for appending to an SFBM (`$add_columns()`).

## bigsnpr 1.9.8

- Fix an issue in `snp_readBGI()` when using an outdated version of package {bit64}.

## bigsnpr 1.9.7

- `snp_cor()` and `bed_cor()` now use less memory.

## bigsnpr 1.9.5

- Remove parameter `info` from `snp_cor()` and `bed_cor()` because this correction is not useful after all.

- `snp_cor()` and `bed_cor()` now return NaNs when e.g. the standard deviation is 0 (and warn about it). Before, these values were not reported (i.e. treated as 0).

## bigsnpr 1.9.4

- You can now return information on all variants with `snp_readBGI()`.

## bigsnpr 1.9.3

- Fix `snp_manhattan()` when non-ordered (chr, pos) are provided.

## bigsnpr 1.9.2

- Enhance function `snp_ancestry_summary()` by allowing to estimate ancestry proportions after PCA projection (instead of directly using the allele frequencies).

## bigsnpr 1.9.1

- Add function `bed_cor()` (similar to `snp_cor()` but with bed files/objects directly).

- Add functions `snp_ld_scores()` and `bed_ld_scores()`.

## bigsnpr 1.9.0

- Add function `snp_ancestry_summary()` to estimate ancestry proportions from a cohort using only its summary allele frequencies.

## bigsnpr 1.8.11

- Add function `snp_scaleAlpha()`, which is similar to `snp_scaleBinom()`, but has a parameter `alpha` that controls the relation between the scaling and the allele frequencies.

## bigsnpr 1.8.10

- Function `snp_cor()` now also uses the upper triangle (`@uplo = "U"`) when the sparse correlation matrix is diagonal, so that it is easier to use with e.g. `as_SFBM()`.

## bigsnpr 1.8.9

- Add parameter `type` in `snp_asGeneticPos()` to also be able to use interpolated genetic maps from [here](https://github.com/joepickrell/1000-genomes-genetic-maps/tree/master/interpolated_from_hapmap).

## bigsnpr 1.8.8

- Add parameter `return_flip_and_rev` to `snp_match()` for whether to return internal boolean variables `"_FLIP_"` and `"_REV_"`.

## bigsnpr 1.8.7

- Add `$perc_kept` in the output of `snp_ldsplit()`, the percentage of initial non-zero values kept within the blocks defined.

## bigsnpr 1.8.6

- Faster `snp_prodBGEN()`.

## bigsnpr 1.8.5

- Add function `snp_prodBGEN()` to compute a matrix product between BGEN files and a matrix (or a vector). This removes the need to read an intermediate FBM object with `snp_readBGEN()` to compute the product. Moreover, when using dosages, they are not rounded to two decimal places anymore.

## bigsnpr 1.8.4

- Trade new parameter `num_iter_change` for a simpler `allow_jump_sign`.

- Change defaults in LDpred2-auto to use 500 burn-in iterations (was 1000 before) followed by 200 iterations (500 before). Such a large number of iterations is usually not really needed.

## bigsnpr 1.8.3 & bigsparser 0.5

- New compact format for SFBMs which should be really useful for LDpred2 (should require about half of memory and be twice as fast). The only thing that you need to change is `as_SFBM(corr0, compact = TRUE)`. Make sure to reinstall {bigsnpr} after updating to {bigsparser} v0.5. 

## bigsnpr 1.8.2

- Prepare for incoming paper on (among other things) improved robustness of LDpred2-auto:
    - add parameter `shrink_corr` to shrink off-diagonal elements of the LD matrix,
    - add parameter `num_iter_change` to control when starting to shrink the variants that change sign too much,
    - also return `corr_est`, the "imputed" correlations between variants and phenotypes, which can be used for post-QCing variants by comparing those to `beta / sqrt(n_eff * beta_se^2 + beta^2)`.

## bigsnpr 1.8.0

- Replace parameter `s` by `delta` in `snp_lassosum2()`. This new parameter `delta` better reflects that the lassosum model also uses L2-regularization (therefore, elastic-net regularization).

- Now detect strong divergence in lassosum2 and LDpred2-grid, and return missing values for the corresponding effect sizes.

## bigsnpr 1.7.4

- Now use a better formula for computing standard errors in `snp_ldsc()` when using blocks with different sizes.

## bigsnpr 1.7.3

- Add parameter `info` to `snp_cor()` to correct correlations when they are computed from imputed dosage data.

## bigsnpr 1.7.2

- Function `snp_readBGEN()` now also returns frequencies and imputation INFO scores.

## bigsnpr 1.7.1

- Add parameter `rsid` to `snp_asGeneticPos()` to also allow matching with rsIDs.

## bigsnpr 1.7.0

- Add function `snp_lassosum2()` to train the lassosum models using the exact same input data as LDpred2.

## bigsnpr 1.6.7

- Add parameter `report_step` in `snp_ldpred2_auto()` to report some of the internal sampling betas.

## bigsnpr 1.6.6

- Fix crash in `snp_readBGEN()` when using BGEN files containing `~`.

## bigsnpr 1.6.5

- Add parameter `thr_r2` in `snp_cor()`.

## bigsnpr 1.6.4

- Remove penalization in `snp_ldsplit()`. Instead, report the best splits for a range of numbers of blocks desired.

## bigsnpr 1.6.3

- Penalization in `snp_ldsplit()` now makes more sense. Also fix a small bug that prevented splitting the last block in some cases.

## bigsnpr 1.6.2

- Add function `snp_ldsplit()` for optimally splitting variants in nearly independent blocks of LD. 

## bigsnpr 1.6.1

- Add option `file.type = "--gzvcf"` for using gzipped VCF in `snp_plinkQC()`.

## bigsnpr 1.6.0

- Finally remove function `snp_assocBGEN()`; prefer reading small parts with `snp_readBGEN()` as a temporary `bigSNP` object and do the association test with e.g. `big_univLinReg()`.

## bigsnpr 1.5.7

- Add function `snp_thr_correct()` for correcting for winner's curse in summary statistics when using p-value thresholding.

## bigsnpr 1.5.6

- Use a better formula for the scale in LDpred2, useful when there are some variants with very large effects (e.g. explaining more than 10% phenotypic variance).

- Simplify LDpred2; there was not really any need for initialization and ordering of the Gibbs sampler.

## bigsnpr 1.5.5

- Add option `return_sampling_betas` in `snp_ldpred2_grid()` to return all sampling betas (after burn-in), which is useful for assessing the uncertainty of the PRS at the individual level (see https://doi.org/10.1101/2020.11.30.403188).

## bigsnpr 1.5.4 & bigsparser 0.4.1

- Faster cross-product with an SFBM, which should make all LDpred2 models faster.

- Also return `$postp_est`, `$h2_init` and `$p_init` in LDpred2-auto.

## bigsnpr 1.5.1

- Add multiple checks in `snp_readBGEN()` to make sure of the expected format.

## bigsnpr 1.5.0

- Add function `snp_fst()` for computing Fst. 

## bigsnpr 1.4.11

- Workaround for error `could not find function "ldpred2_gibbs_auto"`.

## bigsnpr 1.4.9 & bigsparser 0.4.0

- Can now directly do `as_SFBM(corr0)` instead of `bigsparser::as_SFBM(as(corr0, "dgCMatrix"))`. This should also use less memory and be faster.

## bigsnpr 1.4.8

- Add option `sparse` to enable getting also a sparse solution in LDpred2-auto.

## bigsnpr 1.4.7 & bigsparser 0.3.0

- Faster `as_SFBM()`.

- Allow for format `01` or `1` for chromosomes in BGI files.

## bigsnpr 1.4.6

- Fasten `snp_match()`. Also now remove duplicates by default.

## bigsnpr 1.4.3

- Fix a bug when using very large correlation matrices in LDpred2 (although we do not recommend to do so).

## bigsnpr 1.4.2

- All 3 LDpred2 functions now use an SFBM as input format for the correlation matrix.

- Allow for multiple initial values for p in `snp_ldpred2_auto()`.

- Add function `coef_to_liab()` for e.g. converting heritability to the liability scale.

## bigsnpr 1.4.1

- Change default of parameter `alpha` of function `snp_cor()` to `1`.

## bigsnpr 1.4.0

- Add functions `snp_ldpred2_inf()`, `snp_ldpred2_grid()` and `snp_ldpred2_auto()` for running the new LDpred2-inf, LDpred2-grid and LDpred2-auto.

- Add functions `snp_ldsc()` and `snp_ldsc2()` for performing LD score regression.

- Add function `snp_asGeneticPos()` for transforming physical positions to genetic positions.

- Add function `snp_simuPheno()` for simulating phenotypes.

## bigsnpr 1.3.1

- Also use OpenMP for the parallelization of `snp_pcadapt()`, `bed_pcadapt()`, `snp_readBGEN()` and `snp_fastImputeSimple()`.

## bigsnpr 1.3.0

- Parallelization of clumping algorithms has been modified. Before, chromosomes were imputed in parallel. Now, chromosomes are processed sequentially, but computations within each chromosome are performed in parallel thanks to OpenMP. This should prevent major slowdowns for very large samples sizes (due to swapping).

- Use OpenMP to parallelize other functions as well (possibly only sequential until now).

## bigsnpr 1.2.6

- Can now run `snp_cor()` in parallel.

- Parallelization of `snp_fastImpute()` has been modified. Before this version, chromosomes were imputed in parallel. Now, chromosomes are processed sequentially, but computation of correlation between variants and XGBoost models are performed using parallelization.

## bigsnpr 1.2.5

- Add function `snp_subset()` as alias of method `subset()` for subsetting `bigSNP` objects.

## bigsnpr 1.2.4

- Use new class `bed_light` internally to make parallel algorithms faster because they have to transfer less data to clusters. Also define differently functions used in `big_parallelize()` for the same reason.

## bigsnpr 1.2.0

- Use the new implementation of robust OGK Mahalanobis distance in {bigutilsr}.

## bigsnpr 1.1.1

- Fix error `object 'obj.bed' not found` in `snp_readBed2()`.

## bigsnpr 1.1.0

- Cope with new read-only option in {bigstatsr} version >= 1.1.

## bigsnpr 1.0.2

- Add option `backingfile` to `subset.bigSNP()`.

## bigsnpr 1.0.1

- Add option `byrow` to `bed_counts()`.

## bigsnpr 1.0.0

- Add memory-mapping on PLINK (.bed) files with missing values + new functions:
    - `bed()`
    - `bed_MAF()`
    - `bed_autoSVD()`
    - `bed_clumping()`
    - `bed_counts()`
    - `bed_cprodVec()`
    - `bed_pcadapt()`
    - `bed_prodVec()`
    - `bed_projectPCA()`
    - `bed_projectSelfPCA()`
    - `bed_randomSVD()`
    - `bed_scaleBinom()`
    - `bed_tcrossprodSelf()`
    - `download_1000G()`
    - `snp_modifyBuild()`
    - `snp_plinkKINGQC()`
    - `snp_readBed2()`
    - `sub_bed()`
    
- Add 3 parameters to `autoSVD()`: `alpha.tukey`, `min.mac` and `max.iter`.

- Remove option for changing ploidy (that was only partially supported).

- Automatically apply `snp_gc()` to `pcadapt`.

## bigsnpr 0.12.0

- Add `snp_fastImputeSimple()`: fast imputation via mode, mean or sampling according to allele frequencies.

## bigsnpr 0.11.3

- Fix a bug in `snp_readBGEN()` that could not handle duplicated variants or individuals.

## bigsnpr 0.11.1

- When using `snp_grid_PRS()`, it now stores not only the FBM, but also the input parameters as attributes (the whole result basically).

## bigsnpr 0.11.0

- Add 3 SCT functions `snp_grid_*()` to improve from Clumping and Thresholding (preprint coming soon).

- Add `snp_match()` function to match between summary statistics and some SNP information.

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

- **This package won't be on CRAN**. (Okay, it has been back on CRAN since; I was just pissed at BR :D)

## bigsnpr 0.2.4

- No longer download PLINK automatically (because it is a CRAN policy violation).
