/******************************************************************************/

#include <bigstatsr/arma-strict-R-headers.h>
#include <bigstatsr/utils.h>
#include <bigsparser/SFBM.h>

/******************************************************************************/

inline double square(double x) {
  return x * x;
}

inline double gamma_rand(double a, double b) { return ::Rf_rgamma(a, 1 / b); }

inline double mode(double lambda, double omega) {
  double one_minus_lambda = 1 - lambda;
  return omega / (one_minus_lambda + sqrt(square(one_minus_lambda) + square(omega)));
}

double rgig_noshift(double lambda, int check, double omega, double alpha)
{
  double xm,nc,ym,um,s,t,U,V,X;

  t = 0.5*(lambda-1);
  s = 0.25*omega;

  xm = mode(lambda,omega);
  nc = t*log(xm) - s*(xm + 1/xm);
  ym = ((lambda+1) + sqrt((lambda+1)*(lambda+1) + omega*omega))/omega;
  um = exp(0.5*(lambda+1)*log(ym) - s*(ym + 1/ym) - nc);

  do{
    U = um * ::unif_rand();
    V = ::unif_rand();
    X = U/V;
  } while ((log(V)) > (t*log(X) - s*(X+1/X)- nc));

  return (check==1) ? (alpha/X) : (alpha*X);
}

double rgig_shift(double lambda, int check, double omega, double alpha)
{
  double xm,nc,s,t,U,V,X,a,b,c,p,q,fi,fak,y1,y2,uplus,uminus;

  t = 0.5*(lambda-1);
  s = 0.25*omega;

  xm = mode(lambda,omega);
  nc = t*log(xm) - s*(xm + 1/xm);

  a = -(2*(lambda+1)/omega +xm);
  b = (2*(lambda-1)*xm/omega -1);
  c = xm;

  p = b - a*a/3;
  q = (2*a*a*a)/27 - (a*b)/3 + c;

  fi = acos(-q/(2*sqrt(-p*p*p/27)));
  fak = 2*sqrt(-p/3);
  y1 = fak*cos(fi/3) - a/3;
  y2 = fak*cos(fi/3 + (4./3.)*M_PI) - a/3;

  uplus = (y1-xm)*exp(t*log(y1) - s*(y1 + 1/y1) -nc);
  uminus =  (y2-xm)*exp(t*log(y2) - s*(y2 + 1/y2) -nc);

  do {
    U = ::Rf_runif(uminus, uplus);
    V = ::unif_rand();
    X = U / V + xm;
  } while ((X <= 0) || (log(V) > (t * log(X) - s * (X + 1 / X) - nc)));

  return (check == 1) ? (alpha / X) : (alpha * X);
}


double rgig_conc(double lambda, int check, double omega, double alpha)
{
  double A1, A2, A3, Atot,k0,k1,k2,xm,x0,a,U,V,X,hx;

  if (lambda >= 1 || omega > 1)
    Rcpp::stop("Invalid parameters: lambda or omega");

  xm = mode(lambda,omega);
  x0 = omega/(1-lambda);

  k0 = exp((lambda-1)*log(xm) - 0.5*omega*(xm + 1/xm));
  A1 = k0*x0;

  if (x0 >= 2/omega) {
    k1 = 0;
    A2 = 0;
    k2 = pow(x0,lambda-1);
    A3 = k2*2*exp(-omega*x0/2)/omega;
  } else {
    k1 = exp(-omega);
    A2 = (lambda==0) ? (k1*log(2/(omega*omega))) :
      ((k1/lambda)*(pow(2/omega,lambda) - pow(x0,lambda)));
    k2 = pow(2/omega,lambda-1);
    A3 = k2*2*exp(-1.0)/omega;
  }

  Atot = A1 + A2 + A3;

  do {

    V = Atot * ::unif_rand();

    if (V <= A1) {
      X = x0 * V / A1;
      hx = k0;
    } else {

      V -= A1;

      if (V <= A2) {
        if (lambda == 0) {
          X = omega * exp(exp(omega)*V);
          hx = k1 / X;
        } else {
          X = pow(pow(x0, lambda) + (lambda / k1 * V), 1/lambda);
          hx = k1 * pow(X, lambda-1);
        }
      } else {
        V -= A2;
        a = (x0 > 2/omega) ? x0 : 2/omega;
        X = -2/omega * log(exp(-omega/2 * a) - omega/(2*k2) * V);
        hx = k2 * exp(-omega/2 * X);
      }
    }

    U = hx * ::unif_rand();

  } while(log(U) > (lambda-1)*log(X) - omega/2 * (X+1/X));

  return (check==1) ? (alpha/X) : (alpha*X);
}

double rgig(double lambda, double chi, double psi) {

  // code adapted from package {qbld}
  // expecting chi > 0 and psi > 0
  // in my application, lambda is always the same

  int check = 0;

  if (lambda < 0) {
    lambda = -lambda;
    check = 1;
  }

  double alpha = sqrt(psi / chi);
  double omega = sqrt(psi * chi);

  if (lambda > 2 || omega > 3) {
    return rgig_shift(lambda, check, omega, alpha);          // RoU shift
  } else if (lambda >= (1 - 2.25 * chi * psi) || omega > 0.2) {
    return rgig_noshift(lambda, check, omega, alpha);        // RoU no shift
  } else {
    return rgig_conc(lambda, check, omega, alpha);           // log-concave
  }
}

/******************************************************************************/

// [[Rcpp::export]]
arma::vec prscs_gibbs_one(XPtr<SFBM> sfbm,
                          const arma::vec& beta_hat,
                          const NumericVector& n_vec,
                          double a,
                          double b,
                          double phi,
                          double max_psi,
                          int burn_in,
                          int num_iter) {

  double n = Rcpp::max(n_vec);
  int m = beta_hat.size();
  arma::vec curr_beta(m, arma::fill::zeros);
  arma::vec avg_beta(m, arma::fill::zeros);
  arma::vec psi(m, arma::fill::ones);

  double gap0 = arma::dot(beta_hat, beta_hat);

  double cur_h2_est = 0, inv_var_env = 1, a2 = a - 0.5, b2 = a + b;

  for (int k = -burn_in; k < num_iter; k++) {

    double gap = 0;
    double second_part = 0, third_part = 0;

    for (int j = 0; j < m; j++) {

      double dotprod = sfbm->dot_col(j, curr_beta);
      double resid = beta_hat[j] - dotprod;
      gap += resid * resid;
      double res_beta_hat_j = curr_beta[j] + resid;

      double C1 = 1 / (1 + 1 / psi[j]);
      double C2 = res_beta_hat_j * C1;
      if (k > 0) avg_beta[j] += C2;

      double diff = -curr_beta[j];
      double samp_beta = ::Rf_rnorm(C2, sqrt(C1 / n_vec[j] / inv_var_env));
      curr_beta[j] = samp_beta;
      diff += samp_beta;
      cur_h2_est += diff * (2 * dotprod + diff);
      second_part += samp_beta * beta_hat[j];
      third_part += square(samp_beta) / psi[j];
    }

    if (gap > gap0) { avg_beta.fill(NA_REAL); return avg_beta; }

    double other_part = std::max(0.0, 1 - 2 * second_part + cur_h2_est);
    inv_var_env = gamma_rand((n + m) / 2, n / 2 * (third_part + other_part));
    // inv_var_env = 1;

    // Rcout << (k + 1) << " : " << 1 / inv_var_env << " / " << cur_h2_est << std::endl;

    for (int j = 0; j < m; j++) {
      double delta_j = gamma_rand(b2, psi[j] + phi);
      psi[j] = rgig(a2, 2 * delta_j, n_vec[j] * square(curr_beta[j]) * inv_var_env);
      if (psi[j] > max_psi) psi[j] = max_psi;
    }
  }

  return avg_beta / num_iter;
}

/******************************************************************************/

// [[Rcpp::export]]
arma::mat prscs_gibbs(Environment corr,
                      const NumericVector& beta_hat,
                      const NumericVector& n_vec,
                      const NumericVector& a,
                      const NumericVector& b,
                      const NumericVector& phi,
                      double max_psi,
                      int burn_in,
                      int num_iter,
                      int ncores) {

  XPtr<SFBM> sfbm = corr["address"];

  int m = beta_hat.size();
  myassert_size(sfbm->nrow(), m);
  myassert_size(sfbm->ncol(), m);
  myassert_size(n_vec.size(), m);

  int K = a.size();
  myassert_size(b.size(),   K);
  myassert_size(phi.size(), K);

  arma::mat res(m, K);

  #pragma omp parallel for schedule(dynamic, 1) num_threads(ncores)
  for (int k = 0; k < K; k++) {

    // if (k % ncores == 0)
    //   Rcout << "Starting with params " << k + 1 << " / " << K << std::endl;

    arma::vec res_k = prscs_gibbs_one(
      sfbm, beta_hat, n_vec, a[k], b[k], phi[k], max_psi, burn_in, num_iter);

    #pragma omp critical
    res.col(k) = res_k;
  }

  return res;
}

/******************************************************************************/
