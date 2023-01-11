corr0 <- readRDS("~/../Documents/LD_with_blocks_chr22.rds")
ind <- 1:1000 + 3000
corr2 <- as(Matrix::drop0(corr0[ind, ind], tol = sqrt(0.02))^2, "dgCMatrix")

Rcpp::sourceCpp("tmp-tests/test-reorder.cpp")

n <- ncol(corr2)
pos <- 1:n - 1L

for (rep in 1:20) {
  print(rep)
  for (j in 1:ncol(corr2)) {
    if (which_best_place(corr2@p, corr2@i, corr2@x, pos, j - 1L)) {
      cat("Moved", j, "\n")
    }
  }
}

pos <- reorder(corr2, nb_iter = 20)
get_score(corr2@p, corr2@i, corr2@x, pos) # 917809.7

pos
sum(duplicated(pos))
range(pos)

ord <- order(pos)
cowplot::plot_grid(
  Matrix::image(corr2),
  Matrix::image(corr2[ord, ord]),
  nrow = 1
)

Rcpp::sourceCpp("tmp-tests/test-qapSA.cpp")
n <- ncol(corr2)
perm <- qapSA(outer(1:n, 1:n, '-')^2, as.matrix(corr2), 1:n - 1L)
get_score(corr2@p, corr2@i, corr2@x, order(perm)) # 3024962 -> 769203 (smaller init t)
perm2 <- qapSA(outer(1:n, 1:n, '-')^2, as.matrix(corr2), order(pos) - 1L)
get_score(corr2@p, corr2@i, corr2@x, order(perm2)) # 3110880

perm3 <- qap::qap(outer(1:n, 1:n, '-')^2, as.matrix(corr2))
get_score(corr2@p, corr2@i, corr2@x, order(perm3)) # 2694268

Rcpp::sourceCpp("tmp-tests/test-qapSA2.cpp")
perm4 <- qapSA2(outer(1:n, 1:n, '-')^2, as.matrix(corr2), 1:n - 1L, ft = 0.8)
get_score(corr2@p, corr2@i, corr2@x, order(perm4)) # 770309


sum(duplicated(perm))
pos <- order(perm)
(prev_score <- get_score(corr2@p, corr2@i, corr2@x, pos))
ind <- sample(n, 2, replace = FALSE); k1 <- ind[1]; k2 <- ind[2]

delta1 <- sum(corr2[, k1] * ((pos[k2] - pos)^2 - (pos[k1] - pos)^2)) +
  sum(corr2[, k2] * ((pos[k1] - pos)^2 - (pos[k2] - pos)^2)) +
  2 * corr2[k1, k2] * (pos[k1] - pos[k2])^2

vec_diff <- ((k1 - seq_along(perm))^2 - (k2 - seq_along(perm))^2) *
  (corr2[perm + 1L, perm[k1] + 1L] - corr2[perm + 1L, perm[k2] + 1L])
delta2 <- 2 * sum(vec_diff[-ind])

vec_diff2 <- (pos[k1] - pos)^2 - (pos[k2] - pos)^2
delta3 <- sum(((corr2[, k2] - corr2[, k1]) * vec_diff2)[-ind])

delta4 <- sum((corr2[, k2] - corr2[, k1]) * vec_diff2) -
  2 * (1 - corr2[k1, k2]) * (pos[k1] - pos[k2])^2

new_pos <- pos; new_pos[ind] <- pos[rev(ind)]
new_score <- get_score(corr2@p, corr2@i, corr2@x, new_pos)
c(new_score - prev_score, 2 * delta1, delta2, 2 * delta3, 2 * delta4)


Rcpp::sourceCpp("tmp-tests/test-reorder2.cpp")
reorder2(corr2, max_dist = 500, sample_index = TRUE, sample_best = TRUE)
reorder2(corr2, max_dist = 500, sample_index = TRUE, sample_best = FALSE)
reorder2(corr2, max_dist = 500, sample_index = FALSE, sample_best = TRUE)
reorder2(corr2, max_dist = 500, sample_index = FALSE, sample_best = FALSE)

sort(replicate(10, {
  res <- reorder2(corr2, max_dist = 500, sample_index = TRUE, sample_best = TRUE)
  attr(res, "score")
}))
# 766990.0 767421.8 767817.2 767873.7 767937.0 767962.6 768211.8 768561.2 768699.3 770370.3
# 123703115 123770008 123918361 123964874 124742130 125032728 125251736 125298508 125471343 125853259

sort(replicate(10, {
  res <- reorder2(corr2, max_dist = 500, sample_index = FALSE, sample_best = TRUE)
  attr(res, "score")
}))
# 760726.5 765901.1 766268.9 766701.4 766869.0 767952.4 768276.7 768651.6 772485.0 772875.0
# 121290245 121371130 121927334 121948517 121955954 121982183 122086031 122099452 122183395 122784158

sort(replicate(10, {
  res <- reorder2(corr2, max_dist = 500, sample_index = TRUE, sample_best = FALSE)
  attr(res, "score")
}))
# 767099.6 767632.9 767646.3 767855.2 768053.5 768175.2 768405.9 769088.0 769204.3 769984.0
# 123642370 123832146 124933702 124937515 125064954 125074330 125258706 125278039 125883196 125917158

res <- reorder2(corr2, max_dist = 500, sample_index = FALSE, sample_best = FALSE)
attr(res, "score")
# 774128.6
# 122121214


Rcpp::sourceCpp("tmp-tests/test-reorder3.cpp")
pos3 <- reorder3(corr2, max_dist = 50)
