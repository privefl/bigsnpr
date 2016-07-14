#'@title A class for analysis of massive SNP arrays.
#'@description A named list with at least 5 slots:\itemize{
#'\item genotypes: a filebacked \code{\link[bigmemory]{big.matrix}}
#'of type \code{char} representing genotypes.\cr
#'Each element is either 0, 1, 2 or \code{NA}. Rows are individuals and columns are SNPs.
#'\item fam: a \code{data.frame} giving some information on the SNPs.
#'\item map: a \code{data.frame} giving some information on the individuals.
#'\item backingfile: the root name for the backing file(s) for the cache of
#'the object.
#'\item backingpath: the path to the directory containing the file backing cache.
#'}
#'@name bigSNP-class
#'@aliases bigSNP-class bigSNP
#'@keywords class
#'@seealso \code{\link{readplink}}
NULL
