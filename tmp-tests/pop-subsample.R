hist(ldist, ylim = c(0, 1e4))
hist(ldist[na.omit(match(ind_used, ind.row))], add = TRUE,
     col = scales::alpha("red", 0.4))
hist(ldist[ind.sub], add = TRUE, col = scales::alpha("blue", 0.4))

set.seed(4)
ind.sub3 <- sample(length(ldist), 50e3, prob = sqrt(exp(ldist)))
hist(ldist[ind.sub3], add = TRUE, col = scales::alpha("green", 0.4))

system.time(
  obj.svd4 <- bed_autoSVD(obj.bed, ind.row = ind.row[ind.sub3],
                          k = 50, ncores = nb_cores())
) # 3.5H

plot(obj.svd4)

plot_grid(plotlist = lapply(14:25, function(k) {
  plot(obj.svd4, type = "scores", scores = 2 * k - 1:0, coeff = 0.4) +
    aes(color = pop_UKBB[ind.row[ind.sub3]]) +
    theme(legend.position = "none")
}))

set.seed(5)
ind.sub4 <- c(sample(which(ldist < 4), 10e3), which(ldist > 4))
hist(ldist[ind.sub4], add = TRUE, col = scales::alpha("pink", 0.4))

system.time(
  obj.svd5 <- bed_autoSVD(obj.bed, ind.row = ind.row[ind.sub4],
                          k = 50, ncores = nb_cores())
) # 3.5H

plot(obj.svd5)

plot_grid(plotlist = lapply(14:25, function(k) {
  plot(obj.svd5, type = "scores", scores = 2 * k - 1:0, coeff = 0.4) +
    aes(color = pop_UKBB[ind.row[ind.sub4]]) +
    theme(legend.position = "none")
}))

set.seed(6)
ind.sub5 <- sample(length(ldist), 60e3,
                   prob = -pchisq(exp(ldist), 13, lower.tail = FALSE, log.p = TRUE))
hist(ldist)
hist(ldist[ind.sub5], add = TRUE, col = scales::alpha("pink", 0.4))

system.time(
  obj.svd6 <- bed_autoSVD(obj.bed, ind.row = ind.row[ind.sub5],
                          k = 50, ncores = nb_cores())
) # 3.5H

plot(obj.svd6)

plot_grid(plotlist = lapply(14:25, function(k) {
  plot(obj.svd6, type = "scores", scores = 2 * k - 1:0, coeff = 0.4) +
    aes(color = pop_UKBB[ind.row[ind.sub5]]) +
    theme(legend.position = "none")
}))
