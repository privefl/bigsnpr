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
      const NumericVector& curr_beta) {

    // keep only causals + precompute sum_a (which never changes)
    // S always changes, so it proved useless to cache other sums
    nb = nb_causal;
    a = arma::zeros<arma::vec>(nb);
    b = arma::zeros<arma::vec>(nb);
    sum_a = 0;

    int k = 0;
    int m = curr_beta.size();
    for (int j = 0; j < m; j++) {
      double beta_j = curr_beta[j];
      if (beta_j != 0) {
        a[k] = log_var[j];
        b[k] = beta_j * beta_j;
        sum_a += a[k];
        k++;
      }
    }

    if (k != nb) Rcpp::stop("[BUG] with init of causals in MLE.");
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
