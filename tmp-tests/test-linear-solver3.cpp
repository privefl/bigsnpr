// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Eigen;

// [[Rcpp::export]]
VectorXd test_solver_dense(const MatrixXd& A, const VectorXd& b) {
  MINRES<MatrixXd> minres(A);
  return minres.solve(b);
}

// [[Rcpp::export]]
VectorXd test_solver_dense_sparse(const MatrixXd& A) {
  MINRES<MatrixXd> minres(A);
  SparseVector<double> b(A.cols()); b.coeffRef(0) = 1;
  return minres.solve(b);
}

// [[Rcpp::export]]
VectorXd test_ConjugateGradient(const MatrixXd& A) {
  ConjugateGradient<MatrixXd> solver(A);
  SparseVector<double> b(A.cols()); b.coeffRef(0) = 1;
  return solver.solve(b);
}

// [[Rcpp::export]]
VectorXd test_BiCGSTAB(const MatrixXd& A) {
  BiCGSTAB<MatrixXd> solver(A);
  SparseVector<double> b(A.cols()); b.coeffRef(0) = 1;
  return solver.solve(b);
}


// [[Rcpp::export]]
VectorXd test_GMRES(const MatrixXd& A) {
  GMRES<MatrixXd> solver(A);
  SparseVector<double> b(A.cols()); b.coeffRef(0) = 1;
  return solver.solve(b);
}
