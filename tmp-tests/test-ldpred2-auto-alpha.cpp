/******************************************************************************/

#include <bigstatsr/arma-strict-R-headers.h>
#include <bigstatsr/utils.h>
#include <bigsparser/SFBM.h>

/******************************************************************************/

const double EPSILON = 2e-16;
const double MIN_P  = 1e-5;
const double MIN_H2 = 1e-4;

inline double square(double x) {
  return x * x;
}

/******************************************************************************/

double mle_eq(double s,
              const arma::vec& curr_beta,
              const NumericVector& log_var) {
  int m = log_var.size();
  int nb_causal = 0;
  double sum1 = 0, sum2 = 0, sum3 = 0;
  for (int j = 0; j < m; j++) {
    if (curr_beta[j] != 0) {
      double a = log_var[j];
      double b = square(curr_beta[j]) / ::exp(s * a);
      sum1 += b;
      sum2 += a;
      sum3 += a * b;
      nb_causal++;
    }
  }

  return sum1 * sum2 / nb_causal - sum3;
}

double mle_eq2(double s,
               const arma::vec& curr_beta,
               const NumericVector& log_var) {
  int m = log_var.size();
  int nb_causal = 0;
  double sum1 = 0;
  for (int j = 0; j < m; j++) {
    if (curr_beta[j] != 0) {
      sum1 += square(curr_beta[j]) / ::exp(s * log_var[j]);
      nb_causal++;
    }
  }

  return sum1 / nb_causal;
}

// [[Rcpp::export]]
double MLE(const arma::vec& curr_beta,
           const NumericVector& log_var,
           double ax = -0.5,
           double bx = 1.5,
           double tol = 1e-6,
           double maxit = 20) {

  double a = ax, b = bx, c = a;
  double fa = mle_eq(a, curr_beta, log_var), fb = mle_eq(b, curr_beta, log_var), fc = fa;

  if (fa == 0) return a;
  if (fb == 0) return b;

  while (maxit--) {  // Main iteration loop

    double prev_step = b - a;
    double tol_act, p, q, new_step;

    if (fabs(fc) < fabs(fb)) {
      // Swap data for b to be the best approximation
      a = b;    b = c;    c = a;
      fa = fb;  fb = fc;  fc= fa;
    }
    tol_act = 2 * EPSILON * fabs(b) + tol / 2;
    new_step = (c - b) / 2;

    if (fabs(new_step) <= tol_act || fb == 0) return b;

    /* Decide if the interpolation can be tried	*/
    if (fabs(prev_step) >= tol_act && fabs(fa) > fabs(fb)) {
      double t1,cb,t2;
      cb = c-b;
      if( a==c ) {		/* If we have only two distinct	*/
    /* points linear interpolation	*/
    t1 = fb/fa;		/* can only be applied		*/
    p = cb*t1;
    q = 1 - t1;
      }
      else {			/* Quadric inverse interpolation*/

    q = fa/fc;  t1 = fb/fc;	 t2 = fb/fa;
    p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1) );
    q = (q-1) * (t1-1) * (t2-1);
      }
      if( p>(double)0 )		/* p was calculated with the */
    q = -q;			/* opposite sign; make p positive */
    else			/* and assign possible minus to	*/
    p = -p;			/* q				*/

    if( p < (0.75*cb*q-fabs(tol_act*q)/2) && p < fabs(prev_step*q/2) )	/* and isn't too large	*/
    new_step = p/q;			/* it is accepted
     * If p/q is too large then the
     * bisection procedure can
     * reduce [b,c] range to more
     * extent */
    }

    if( fabs(new_step) < tol_act) {	/* Adjust the step to be not less*/
    if( new_step > (double)0 )	/* than tolerance		*/
    new_step = tol_act;
    else
      new_step = -tol_act;
    }
    a = b;	fa = fb;			/* Save the previous approx. */
    b += new_step;	fb = mle_eq(b, curr_beta, log_var);	/* Do step to a new approxim. */
    if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0) ) {
      /* Adjust c for it to have a sign opposite to that of b */
      c = a;  fc = fa;
    }
  }

  return (fabs(fb) < fabs(fa)) ? b : a;  // failed, returns best solution for now
}

/******************************************************************************/

// [[Rcpp::export]]
List ldpred2_gibbs_auto_alpha(Environment corr,
                              const NumericVector& beta_hat,
                              const NumericVector& beta_init,
                              const IntegerVector& order,
                              const NumericVector& n_vec,
                              const NumericVector& log_var,
                              double p_init,
                              double h2_init,
                              int burn_in,
                              int num_iter,
                              int report_step,
                              bool allow_jump_sign,
                              double shrink_corr,
                              bool verbose = false) {

  XPtr<SFBM> sfbm = corr["address"];

  int m = beta_hat.size();
  myassert_size(sfbm->nrow(), m);
  myassert_size(sfbm->ncol(), m);
  myassert_size(order.size(), m);
  myassert_size(beta_init.size(), m);
  myassert_size(n_vec.size(), m);

  arma::vec curr_beta(beta_init.begin(), m);
  arma::vec avg_beta(m, arma::fill::zeros), avg_postp(m, arma::fill::zeros);
  arma::vec avg_beta_hat(m, arma::fill::zeros);

  arma::mat sample_beta(m, num_iter / report_step, arma::fill::zeros);
  int ind_report = 0, next_k_reported = burn_in - 1 + report_step;

  int num_iter_tot = burn_in + num_iter;
  std::vector<double> p_est(num_iter_tot), h2_est(num_iter_tot), alpha(num_iter_tot);

  double cur_h2_est = arma::dot(curr_beta, sfbm->prod(curr_beta));
  double p = p_init, h2 = h2_init, s_plus_one = 0;
  double avg_p = 0, avg_h2 = 0, sigma2 = h2 / (m * p);

  for (int k = 0; k < num_iter_tot; k++) {

    int nb_causal = 0;
    double inv_odd_p = (1 - p) / p;
    // double sum_scale_freq = 0, sum_xx = 0, sum_x = 0, sum_y = 0, sum_xy = 0;

    for (const int& j : order) {

      double dotprod = sfbm->dot_col(j, curr_beta);
      double res_beta_hat_j = beta_hat[j] + shrink_corr * (curr_beta[j] - dotprod);

      double scale_freq = ::exp(s_plus_one * log_var[j]);
      double C1 = scale_freq * sigma2 * n_vec[j];
      double C2 = 1 / (1 + 1 / C1);
      double C3 = C2 * res_beta_hat_j;
      double C4 = ::sqrt(C2 / n_vec[j]);

      double postp = 1 /
        (1 + inv_odd_p * ::sqrt(1 + C1) * ::exp(-square(C3 / C4) / 2));

      double prev_beta = curr_beta[j];
      double dotprod_shrunk = shrink_corr * dotprod + (1 - shrink_corr) * prev_beta;

      if (k >= burn_in) {
        avg_postp[j]    += postp;
        avg_beta[j]     += C3 * postp;
        avg_beta_hat[j] += dotprod_shrunk;
      }

      double diff = -prev_beta;
      if (postp > ::unif_rand()) {

        double samp_beta = ::Rf_rnorm(C3, C4);

        if (!allow_jump_sign && (samp_beta * prev_beta) < 0) {
          curr_beta[j] = 0;
          if (k >= burn_in) {
            avg_postp[j] -= postp;
            avg_beta[j]  -= C3 * postp;
          }
        } else {
          curr_beta[j] = samp_beta;
          diff += samp_beta;
          nb_causal++;

          // double x_j = log_var[j];
          // sum_x += x_j;
          // sum_xx += square(x_j);
          // double y_j = log(square(samp_beta));
          // sum_y += y_j;
          // sum_xy += y_j * x_j;
          //
          // sum_scale_freq += scale_freq;
        }

      } else {
        curr_beta[j] = 0;
      }

      cur_h2_est += diff * (2 * dotprod_shrunk + diff);
    }

    p = std::max(::Rf_rbeta(1 + nb_causal, 1 + m - nb_causal), MIN_P);
    h2 = std::max(cur_h2_est, MIN_H2);
    // s_plus_one = (sum_xy - sum_x * sum_y / nb_causal) /
    //   (sum_xx - square(sum_x) / nb_causal);
    s_plus_one = std::min(std::max(-0.5, MLE(curr_beta, log_var)), 1.5);
    sigma2 = mle_eq2(s_plus_one, curr_beta, log_var);

    if (verbose) Rcout << k + 1 << ": " << p << " // " << h2 <<
      " // " << (s_plus_one - 1) << std::endl;

    if (k >= burn_in) {
      avg_p  += p;
      avg_h2 += h2;
      if (k == next_k_reported) {
        sample_beta.col(ind_report++) = curr_beta;
        next_k_reported += report_step;
      }
    }
    p_est[k]  = p;
    h2_est[k] = h2;
    alpha[k] = s_plus_one - 1;
  }

  double est_p  = avg_p  / num_iter;
  double est_h2 = avg_h2 / num_iter;
  if (verbose) Rcout << "Overall: " << est_p << " // " << est_h2 << std::endl;

  return List::create(
    _["beta_est"]    = avg_beta     / num_iter,
    _["postp_est"]   = avg_postp    / num_iter,
    _["corr_est"]    = avg_beta_hat / num_iter,
    _["sample_beta"] = sample_beta,
    _["p_est"]       = est_p,
    _["h2_est"]      = est_h2,
    _["path_p_est"]  = p_est,
    _["path_h2_est"] = h2_est,
    _["path_alpha_est"] = alpha);
}

/******************************************************************************/
