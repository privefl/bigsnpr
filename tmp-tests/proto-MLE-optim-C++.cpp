#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <roptim.h>
// [[Rcpp::depends(roptim)]]

using namespace roptim;

class MLE : public Functor {
public:
  MLE(const arma::vec& a_, const arma::vec& b_) : a(a_), b(b_) {
    // precomputing sum_a (which never changes)
    // S always changes, so it seems useless to cache other sums
    sum_a = 0;
    m = a.size();
    for (int j = 0; j < m; j++) sum_a += a[j];
  }

  double operator()(const arma::vec& par) override {

    double S = par[0] + 1;
    double sigma2 = par[1];

    double sum_c = 0;
    for (int j = 0; j < m; j++)
      sum_c += b[j] * ::exp(-S * a[j]);

    return S * sum_a + m * ::log(sigma2) + sum_c / sigma2;
  }

  void Gradient(const arma::vec& par, arma::vec& gr) override {

    double S = par[0] + 1;
    double sigma2 = par[1];

    double sum_c = 0, sum_ac = 0;
    for (int j = 0; j < m; j++) {
      double c_j = b[j] * ::exp(-S * a[j]);
      sum_c  += c_j;
      sum_ac += a[j] * c_j;
    }

    gr = arma::zeros<arma::vec>(2);
    gr[0] = sum_a - sum_ac / sigma2;
    gr[1] = (m - sum_c / sigma2) / sigma2;
  }

private:
  arma::vec a;
  arma::vec b;
  int m;
  double sum_a;
};


// [[Rcpp::export]]
arma::vec test_MLE(const arma::vec& a_, const arma::vec& b_,
                   arma::vec par, const arma::vec& lower, const arma::vec& upper) {

  MLE mle(a_, b_);
  Roptim<MLE> opt("L-BFGS-B");
  opt.set_lower(lower);
  opt.set_upper(upper);
  opt.set_hessian(false);
  opt.control.trace = 0;

  opt.minimize(mle, par);

  Rcpp::Rcout << "-------------------------" << std::endl;
  opt.print();

  return par;
}
