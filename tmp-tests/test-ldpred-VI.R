library(bigsnpr)

chr22 <- snp_attach("../Dubois2010_data/celiac_chr22.rds")
G <- chr22$genotypes$copy(code = c(0, 1, 2, 0, rep(NA, 252)))
dim(G) #  11402 x 4945
big_counts(G, ind.col = 1:10)
CHR <- chr22$map$chromosome
POS <- chr22$map$physical.pos
stats <- big_scale()(G)

POS2 <- snp_asGeneticPos(CHR, POS, dir = "tmp-data/")
plot(POS, POS2, pch = 20)

corr <- runonce::save_run(
  snp_cor(chr22$genotypes, infos.pos = POS2, size = 3 / 1000, ncores = 6),
  file = "tmp-data/corr_chr22.rds"
)

# hist(log10(abs(corr@x)))

# Simu phenotype
set.seed(1)
(M <- round(ncol(G) * 10^-runif(1, 1, 2)))  # 268
simu <- snp_simuPheno(G, h2 = 0.1, M = M)

y <- simu$pheno

# GWAS
set.seed(1)
ind.gwas <- sample(nrow(G), 8e3)
ind.val <- setdiff(rows_along(G), ind.gwas)
gwas <- big_univLinReg(G, y[ind.gwas], ind.train = ind.gwas)
plot(gwas)
plot(gwas, type = "Manhattan")

df_beta <- data.frame(beta = gwas$estim, beta_se = gwas$std.err,
                      n_eff = length(ind.gwas))



# input parameters
n_vec <- df_beta$n_eff
beta_hat <- with(df_beta, beta / sqrt(n_eff * beta_se^2 + beta^2))

mean_ld <- mean(Matrix::colSums(corr^2))

burn_in <- 100
num_iter <- 100

p_init <- 0.05
h2_init <- 0.1


#### LDpred2-auto ####

corr2 <- as.matrix(corr)

m <- length(beta_hat)

{
  curr_beta <- dotprods <- avg_beta <- avg_postp <- avg_beta_hat <-
    rep(0, m)
  num_iter_tot <- burn_in + num_iter
  p_est <- h2_est <- rep(NA_real_, num_iter_tot)

  id <- matrix(0, m, 1)

  cur_h2_est <- 0
  p <- p_init
  h2 <- h2_init
  gap0 <- crossprod(beta_hat)
  shrink_corr <- 0.999

  for (k in seq_len(num_iter_tot)) {

    inv_odd_p = (1 - p) / p
    sigma2 = h2 / (m * p)
    gap = 0
    nb_causal <- 0

    for (j in seq_len(m)) {

      # print(j)

      dotprod = dotprods[j];
      resid = beta_hat[j] - dotprod;
      gap =  gap + resid * resid;
      res_beta_hat_j = beta_hat[j] + shrink_corr * (curr_beta[j] - dotprod);

      C1 = sigma2 * n_vec[j];
      C2 = 1 / (1 + 1 / C1);
      C3 = C2 * res_beta_hat_j;
      C4 = C2 / n_vec[j];

      postp = 1 / (1 + inv_odd_p * sqrt(1 + C1) * exp(-C3 * C3 / C4 / 2));

      dotprod_shrunk = shrink_corr * dotprod + (1 - shrink_corr) * curr_beta[j];

      if (k > burn_in) {
        avg_postp[j]    = avg_postp[j] + postp;
        avg_beta[j]     = avg_beta[j] + C3 * postp;
        avg_beta_hat[j] = avg_beta_hat[j] + dotprod_shrunk;
      }

      if (postp > runif(1)) {

        samp_beta = rnorm(1, C3, sqrt(C4));

        diff = samp_beta - curr_beta[j];
        curr_beta[j] = samp_beta;
        nb_causal <- nb_causal + 1

      } else {
        diff = -curr_beta[j];
        curr_beta[j] = 0;
      }

      if (diff != 0) {
        cur_h2_est = cur_h2_est + diff * (2 * dotprod_shrunk + diff);
        dotprods = dotprods + corr2[, j] * diff
        # id[[j]] <- diff
        # dotprods = dotprods + Matrix::solve(prec2, id)[, 1]
        # id[[j]] <- 0
      }
    }

    if (gap > gap0) stop("Divergence!")

    # p = rbeta(1, 1 + nb_causal / mean_ld, 1 + (m - nb_causal) / mean_ld)
    # h2 = cur_h2_est

    print(c(k, p, h2))
    p_est[k]  = p;
    h2_est[k] = h2;
  }
}

plot(p_est, log = "y", pch = 20); abline(h = M / m, col = "red")
plot(h2_est, pch = 20); abline(h = 0.2, col = "red")

pred <- big_prodVec(G, avg_beta, ind.row = ind.val)
cor(pred, y[ind.val])^2 # 6.1%
# using the normal corr leads to divergence


## Variational inference (based on https://bit.ly/ldpred_VI)

{
  curr_beta <- dotprods <- avg_beta <- avg_postp <- avg_beta_hat <- postp <-
    rep(0, m)
  num_iter_tot <- burn_in + num_iter
  p_est <- h2_est <- rep(NA_real_, num_iter_tot)

  id <- matrix(0, m, 1)

  cur_h2_est <- 0
  p <- p_init
  h2 <- h2_init
  gap0 <- crossprod(beta_hat)
  shrink_corr <- 0.999

  tol <- median(abs(beta_hat)) * 1e-2

  for (k in seq_len(num_iter_tot)) {

    print(k)
    stop <- TRUE

    inv_odd_p = (1 - p) / p
    sigma2 = h2 / (m * p)
    gap = 0
    nb_causal <- 0

    prev_beta <- curr_beta

    for (j in seq_len(m)) {

      # print(j)

      dotprod = dotprods[j];
      resid = beta_hat[j] - dotprod;
      gap =  gap + resid * resid;
      res_beta_hat_j = beta_hat[j] + shrink_corr * (curr_beta[j] - dotprod);

      C1 = sigma2 * n_vec[j];
      C2 = 1 / (1 + 1 / C1);
      C3 = C2 * res_beta_hat_j;
      C4 = C2 / n_vec[j];

      postp[j] = 1 / (1 + inv_odd_p * sqrt(1 + C1) * exp(-C3 * C3 / C4 / 2));

      dotprod_shrunk = shrink_corr * dotprod + (1 - shrink_corr) * curr_beta[j];

      if (k > burn_in) {
        avg_postp[j]    = avg_postp[j] + postp[j];
        avg_beta[j]     = avg_beta[j] + C3 * postp[j];
        avg_beta_hat[j] = avg_beta_hat[j] + dotprod_shrunk;
      }

      new_beta <- C3 * postp[j]
      diff <- new_beta - prev_beta[j]
      if (abs(diff) > tol) {
        curr_beta[j] = new_beta
        dotprods = dotprods + corr2[, j] * diff
        stop <- FALSE
      }
    }

    if (stop) break
    # nb_causal <- sum(postp)
    # p = rbeta(1, 1 + nb_causal / mean_ld, 1 + (m - nb_causal) / mean_ld)
  }
}

pred3 <- big_prodVec(G, curr_beta, ind.row = ind.val)
cor(pred3, y[ind.val])^2 # 5.8%
