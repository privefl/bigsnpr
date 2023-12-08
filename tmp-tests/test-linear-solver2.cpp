// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

class MatrixReplacement;
using Eigen::SparseMatrix;

namespace Eigen {
namespace internal {
// MatrixReplacement looks-like a SparseMatrix, so let's inherits its traits:
template<>
struct traits<MatrixReplacement> :  public Eigen::internal::traits<Eigen::SparseMatrix<double> >
{};
}
}

// Example of a matrix-free wrapper from a user type to Eigen's compatible type
// For the sake of simplicity, this example simply wrap a Eigen::SparseMatrix.
class MatrixReplacement : public Eigen::EigenBase<MatrixReplacement> {
public:
  // Required typedefs, constants, and method:
  typedef double Scalar;
  typedef double RealScalar;
  typedef int StorageIndex;
  enum {
    ColsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic,
    IsRowMajor = false
  };

  Index rows() const { return mp_mat->rows(); }
  Index cols() const { return mp_mat->cols(); }

  template<typename Rhs>
  Eigen::Product<MatrixReplacement,Rhs,Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
    return Eigen::Product<MatrixReplacement,Rhs,Eigen::AliasFreeProduct>(*this, x.derived());
  }

  // Custom API:
  MatrixReplacement() : mp_mat(0) {}

  void attachMyMatrix(const SparseMatrix<double> &mat) {
    mp_mat = &mat;
  }
  const SparseMatrix<double> my_matrix() const { return *mp_mat; }

private:
  const SparseMatrix<double> *mp_mat;
};


// Implementation of MatrixReplacement * Eigen::DenseVector though a specialization of internal::generic_product_impl:
namespace Eigen {
namespace internal {

template<typename Rhs>
struct generic_product_impl<MatrixReplacement, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
  : generic_product_impl_base<MatrixReplacement,Rhs,generic_product_impl<MatrixReplacement,Rhs> >
{
  typedef typename Product<MatrixReplacement,Rhs>::Scalar Scalar;

  template<typename Dest>
  static void scaleAndAddTo(Dest& dst, const MatrixReplacement& lhs, const Rhs& rhs, const Scalar& alpha)
  {
    // This method should implement "dst += alpha * lhs * rhs" inplace,
    // however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
    assert(alpha==Scalar(1) && "scaling is not implemented");
    EIGEN_ONLY_USED_FOR_DEBUG(alpha);

    dst.noalias() += lhs.my_matrix() * rhs;
  }
};

}
}

// [[Rcpp::export]]
Eigen::VectorXd test_solver(Eigen::SparseMatrix<double>& S, Eigen::VectorXd& b) {

  MatrixReplacement A;
  A.attachMyMatrix(S);
  Eigen::MINRES<MatrixReplacement, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner> minres;
  minres.compute(A);
  Eigen::VectorXd x = minres.solve(b);

  return x;
}
