#ifndef OPTIM_MLE_ALPHA_H
#define OPTIM_MLE_ALPHA_H

/******************************************************************************/

#include <bigstatsr/arma-strict-R-headers.h>
#include <roptim.h>

using namespace roptim;

/******************************************************************************/

class MLE : public Functor {
public:
  // constructor
  MLE(int nb_causal,
      const NumericVector& log_var,
      const NumericVector& curr_beta,
      bool boot = false) {

    // keep only causals + precompute sum_a (which never changes)
    // S always changes, so it proved useless to cache other sums
    nb = nb_causal;
    a = arma::zeros<arma::vec>(nb);
    b = arma::zeros<arma::vec>(nb);

    int k = 0;
    int m = curr_beta.size();
    for (int j = 0; j < m; j++) {
      double beta_j = curr_beta[j];
      if (beta_j != 0) {
        a[k] = log_var[j];
        b[k] = beta_j * beta_j;
        k++;
      }
    }
    if (k != nb) Rcpp::stop("[BUG] with init of causals in MLE.");

    if (boot) {  // useful to take uncertainty into account
      arma::uvec ind(nb);
      for (int k = 0; k < nb; k++)
        ind[k] = nb * unif_rand();

      a = a(ind);
      b = b(ind);
    }

    sum_a = arma::sum(a);
  }

  // objective function
  double operator()(const arma::vec& par) override {

    double S = par[0];
    double sigma2 = par[1];

    double sum_c = 0;
    for (int k = 0; k < nb; k++)
      sum_c += b[k] * ::exp(-S * a[k]);

    return S * sum_a + nb * ::log(sigma2) + sum_c / sigma2;
  }

  // partial derivatives
  void Gradient(const arma::vec& par, arma::vec& gr) override {

    double S = par[0];
    double sigma2 = par[1];

    double sum_c = 0, sum_ac = 0;
    for (int k = 0; k < nb; k++) {
      double c_k = b[k] * ::exp(-S * a[k]);
      sum_c  += c_k;
      sum_ac += a[k] * c_k;
    }

    gr = arma::zeros<arma::vec>(2);
    gr[0] = sum_a - sum_ac / sigma2;
    gr[1] = (nb - sum_c / sigma2) / sigma2;
  }

private:
  arma::vec a;
  arma::vec b;
  int nb;
  double sum_a;
};

/******************************************************************************/

#endif // OPTIM_MLE_ALPHA_H
