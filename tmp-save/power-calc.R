# I am trying to create graphs that show power to detect variants of different effect sizes (OR) for MAFs of 1, 5, 10, 20% for common variants, and 0.1% for rare variants, all at credible alpha values. I have 7000 controls and 700 cases.
prevalence <- 0.1
alpha <- 5e-8
num.cases <- 700
num.controls <- 7000
num.total <- num.cases + num.controls
prop.cases <- num.cases / num.total

grid <- expand.grid(
  MAF = c(0.1, 1, 5, 10, 20) / 100,
  OR = seq(1.1, 1.5, by = 0.1)
)

grid$power <- purrr::pmap_dbl(grid, function(MAF, OR) {
  gap::pbsize2(
    N = num.total,
    fc = prop.cases,
    alpha = alpha,
    gamma = OR,
    p = MAF,
    kp = prevalence,
    model = "additive"
  )
})

library(ggplot2)
ggplot(grid, aes(OR, power)) +
  bigstatsr::theme_bigstatsr() +
  geom_line() +
  facet_wrap(~ MAF)

