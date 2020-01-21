################################################################################

context("MAX3")

opt.save <- options(bigstatsr.typecast.warning = FALSE)

# from p173 of Zheng2012 -> two errors on p1
tab <- c( 25, 283, 864,  10, 218, 929, 0.0150, 0.2151,
         223, 598, 351, 301, 579, 277, 0.2250, 0.5054,
          27, 283, 861,  11, 206, 939, 0.0163, 0.2101,
          10, 180, 955,  14, 272, 854, 0.0105, 0.1978,
          50, 477, 608,  99, 408, 628, 0.0656, 0.3899,
          18, 316, 777,  26, 220, 862, 0.0198, 0.2416,
         250, 543, 352, 170, 538, 433, 0.1837, 0.4729,
         187, 605, 353, 249, 496, 396, 0.1907, 0.4816,
         242, 546, 357, 165, 537, 440, 0.1780, 0.4735)
dim(tab) <- c(8, 9)
p0 <- (tab[1, ] + tab[4, ]) / colSums(tab[1:6, ])
p1 <- (tab[2, ] + tab[5, ]) / colSums(tab[1:6, ])

test_that("Same pi as reported", {
  expect_equal(round(p0, 4), tab[7, ])
  expect_equal(round(p1, 4), tab[8, ])
})

res <- c(4.080, 4.468, 4.694, 4.999, 4.153,
         4.214, 4.773, 3.341, 4.759) # p174

################################################################################

N.cases <- max(colSums(tab[1:3, ]))
N.casesNA <- N.cases - colSums(tab[1:3, ])
N.controls <- max(colSums(tab[4:6, ]))
N.controlsNA <- N.controls - colSums(tab[4:6, ])

fake <- snp_fake(N.cases + N.controls, 9)
fake$fam$affection <- c(rep(2, N.cases), rep(1, N.controls))

G <- fake$genotypes
for (j in 1:9) {
  cases <- c(rep(0, tab[1, j]),
             rep(1, tab[2, j]),
             rep(2, tab[3, j]),
             rep(NA, N.casesNA[j]))
  controls <- c(rep(0, tab[4, j]),
                rep(1, tab[5, j]),
                rep(2, tab[6, j]),
                rep(NA, N.controlsNA[j]))
  vals <- c(sample(cases), sample(controls))
  G[, j] <- as.raw(replace(vals, is.na(vals), 3))
}

################################################################################

test <- snp_MAX3(Gna = fake$genotypes, y01.train = fake$fam$affection - 1)

test_that("Same values for statistics", {
  expect_equal(round(sqrt(test$score), 3), res)
})

test_that("Same ranks", {
  expect_equal(rank(test$score), rank(-predict(test)))
})

################################################################################

options(opt.save)

################################################################################
