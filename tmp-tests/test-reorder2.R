library(magrittr)  # for the %>%
source_cpp_url <- function(url) {
    runonce::download_file(url, dir = "tmp-tests") %>%
    Rcpp::sourceCpp()
}

corr0 <- "https://figshare.com/ndownloader/files/24928052" %>%
  runonce::download_file(dir = "tmp-data") %>%
  readRDS()

ind <- 1:500 + 2000
corr0_sq <- corr0[ind, ind] %>%
  Matrix::drop0(tol = 0.1) %>%
  .^2 %>%
  as("generalMatrix")  # not symmetric storing anymore; easier to work with

Matrix::image(corr0_sq)  # there should be a grey gradient depending on r2


source_cpp_url("https://raw.githubusercontent.com/privefl/bigsnpr/master/tmp-tests/test-reorder2.cpp")
source_cpp_url("https://raw.githubusercontent.com/privefl/bigsnpr/master/tmp-tests/test-qapSA.cpp")
n <- ncol(corr0_sq)
get_score(corr0_sq@p, corr0_sq@i, corr0_sq@x, 1:n)         # 8682899
perm <- qapSA(outer(1:n, 1:n, '-')^2, as.matrix(corr0_sq), 1:n - 1L) + 1L
get_score(corr0_sq@p, corr0_sq@i, corr0_sq@x, order(perm)) # 6492697

perm3 <- qap::qap(outer(1:n, 1:n, '-')^2, as.matrix(corr0_sq))
attr(perm3, "obj")  # 10212726 -> do not remember why so large
get_score(corr0_sq@p, corr0_sq@i, corr0_sq@x, order(perm3)) # 10212726


source_cpp_url("https://raw.githubusercontent.com/privefl/bigsnpr/master/tmp-tests/test-reorder2.cpp")
# tries to make some permutation i <--> j
reorder2(corr0_sq, max_dist = 500, sample_index = TRUE, sample_best = TRUE)
reorder2(corr0_sq, max_dist = 500, sample_index = TRUE, sample_best = FALSE)
reorder2(corr0_sq, max_dist = 500, sample_index = FALSE, sample_best = TRUE)
reorder2(corr0_sq, max_dist = 500, sample_index = FALSE, sample_best = FALSE)
# all the sampling does not seem very useful + not deterministic

pos2 <- reorder2(corr0_sq, max_dist = 500, sample_index = TRUE, sample_best = TRUE)


source_cpp_url("https://raw.githubusercontent.com/privefl/bigsnpr/master/tmp-tests/test-reorder3.cpp")
# tries to move i to some other position
pos3 <- reorder3(corr0_sq, max_dist = 50)
get_score(corr0_sq@p, corr0_sq@i, corr0_sq@x, pos3) # 6540525

ord2 <- order(pos2)
ord3 <- order(pos3)
library(latticeExtra)
c(before = Matrix::image(corr0_sq),
  after_qapSA = Matrix::image(corr0_sq[perm, perm]),
  after_reorder2 = Matrix::image(corr0_sq[ord2, ord2]),
  layout = c(3, 1))
