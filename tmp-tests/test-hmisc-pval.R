# more in corr than in true$P
# so p-values are lower than it should
# N is bigger than it should?

n <- true$n[50, 49]
r <- true$r[50, 49]
r <- 0.2
t <- r * sqrt((n-2)/(1-r^2))
2 * pt(abs(t), n-2, lower.tail = FALSE)
true$P[50, 49]

2 * p(abs(t)) = 0.05
p(abs(t)) = 0.05 / 2
abs(t) = q(0.05 / 2)
