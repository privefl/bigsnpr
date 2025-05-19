N <- 10000
# if we assume an AA0 distribution similar to e.g. breast cancer
curve(dnorm(x, mean = 40, sd = 10), from = 0, to = 100)
# if we have some full liability
y <- rnorm(N)
# and a lifetime prevalence of say 15%, the threshold t for being a case is:
K <- 0.15
t <- pnorm(K, lower.tail = FALSE)
y01 <- (y > t)
# and the age of onset based on the liability is
hist(p_y <- ifelse(y < t, NA, truncnorm::ptruncnorm(y, t)))
hist(aao <- qnorm(p_y, mean = 40, sd = 10, lower.tail = FALSE))
plot(y, aao)
# you can also try to add some noise to the AAO
aao_noise <- aao + runif(N, max = 10)
plot(y, aao_noise)

hazard <- exp(-y)  # Transform liability into a hazard ratio
hist(aao2 <- rweibull(N, shape = 1, scale = hazard))
library(ggplot2)
qplot(aao, aao2) +
  # coord_cartesian(ylim = c(NA, 5)) +
  scale_y_log10()
