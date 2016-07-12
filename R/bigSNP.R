#'@title R Package for analysis of massive SNP arrays.
#'@description A named list with at least 3 slots:\itemize{
#'\item genotypes: a filebacked \code{\link[bigmemory]{big.matrix}}
#'of type \code{char} representing genotypes.\cr
#'Each element is either 0, 1, 2 or NA. Rows are individuals and columns are SNPs.
#'\item fam: a \code{data.frame} giving some information on the SNPs.
#'\item map: a \code{data.frame} giving some information on the individuals.
#'}
#'Others slots will be added by methods of \code{bigSNP},
#'for instance \code{\link[bigsnpr]{impute}}.
#'@name bigSNP-class
#'@aliases bigSNP-class bigSNP
#'@keywords class
#'@seealso \code{\link{readplink}}
NULL
