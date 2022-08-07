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
  MLE(const std::vector<int>& ind_causal,
      const NumericVector& log_var,
      const NumericVector& curr_beta,
      bool boot = false) {

    // keep only causals + precompute sum_a (which never changes)
    // alpha always changes, so it proved useless to cache other sums
    nb = ind_causal.size();
    a = arma::zeros<arma::vec>(nb);
    b = arma::zeros<arma::vec>(nb);

    for (int k = 0; k < nb; k++) {
      int k2 = boot ? nb * unif_rand() : k;
      int j = ind_causal[k2];
      a[k] = log_var[j];
      b[k] = curr_beta[j] * curr_beta[j];
    }

    sum_a = arma::sum(a);
  }

  // objective function
  double operator()(const arma::vec& par) override {

    double alpha_plus_one = par[0];
    double sigma2 = par[1];

    double sum_c = 0;
    for (int k = 0; k < nb; k++)
      sum_c += b[k] * ::exp(-alpha_plus_one * a[k]);

    return alpha_plus_one * sum_a + nb * ::log(sigma2) + sum_c / sigma2;
  }

  // partial derivatives
  void Gradient(const arma::vec& par, arma::vec& gr) override {

    double alpha_plus_one = par[0];
    double sigma2 = par[1];

    double sum_c = 0, sum_ac = 0;
    for (int k = 0; k < nb; k++) {
      double c_k = b[k] * ::exp(-alpha_plus_one * a[k]);
      sum_c  += c_k;
      sum_ac += a[k] * c_k;
    }

    gr[0] = sum_a - sum_ac / sigma2;
    gr[1] = (nb - sum_c / sigma2) / sigma2;
  }

private:
  int nb;
  arma::vec a;
  arma::vec b;
  double sum_a;
};

/******************************************************************************/

#endif // OPTIM_MLE_ALPHA_H
