/******************************************************************************/

// [[Rcpp::depends(RcppEigen)]]
#include <bigsparser/EigenMatrixReplacement.h>
typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

/******************************************************************************/

// [[Rcpp::export]]
void sp_solve_sym_eigen(Rcpp::Environment X,
                        const Eigen::VectorXd& b,
                        const Eigen::VectorXd& add_to_diag,
                        double tol,
                        int maxiter) {

  Rcpp::XPtr<SFBM> sfbm = X["address"];
  MatrixReplacement A(sfbm, add_to_diag);

  // Solve Ax = b using iterative solvers with matrix-free version
  // Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower | Eigen::Upper,
  //                          Eigen::IdentityPreconditioner> solver;
  Eigen::MINRES<MatrixReplacement, Eigen::Lower | Eigen::Upper,
                Eigen::IdentityPreconditioner> solver;

  solver.setTolerance(tol);
  solver.setMaxIterations(maxiter);

  solver.compute(A);
  // Eigen::VectorXd x = solver.solve(b);

  double eps = solver.error();
  // std::cout << eps << " in " << solver.iterations() << " iterations." << std::endl;

  if (std::isnan(eps))
    Rcpp::stop("Solver failed.");
  if (eps > tol)
    Rcpp::warning("Estimated error: %s.", eps);
}

/******************************************************************************/

// [[Rcpp::export]]
Rcpp::NumericVector sp_solve_sym_eigen2(const SpMat& A,
                                        const Eigen::VectorXd& b,
                                        double tol,
                                        int maxiter) {

  Eigen::MINRES<SpMat> solver;

  solver.setTolerance(tol);
  solver.setMaxIterations(maxiter);

  solver.compute(A);
  Eigen::VectorXd x = solver.solve(b);

  double eps = solver.error();
  // std::cout << eps << " in " << solver.iterations() << " iterations." << std::endl;

  if (std::isnan(eps))
    Rcpp::stop("Solver failed.");
  if (eps > tol)
    Rcpp::warning("Estimated error: %s.", eps);

  return Rcpp::wrap(x);
}

/******************************************************************************/

// [[Rcpp::export]]
Rcpp::NumericVector sp_solve_sym_eigen3(const SpMat& A,
                                        const Eigen::VectorXd& b,
                                        double tol,
                                        int maxiter) {

  Eigen::MINRES<SpMat > solver;

  solver.setTolerance(tol);
  solver.setMaxIterations(maxiter);

  solver.compute(A);
  Eigen::SparseVector<double> sp_b = b.sparseView(); // create sparse vector
  Eigen::VectorXd x = solver.solve(sp_b);

  double eps = solver.error();
  // std::cout << eps << " in " << solver.iterations() << " iterations." << std::endl;

  if (std::isnan(eps))
    Rcpp::stop("Solver failed.");
  if (eps > tol)
    Rcpp::warning("Estimated error: %s.", eps);

  return Rcpp::wrap(x);
}

/******************************************************************************/

// [[Rcpp::export]]
Rcpp::NumericVector sp_solve_sym_eigen4(const SpMat& A,
                                        const Eigen::VectorXd& b,
                                        int nb_solve = 1) {

  Eigen::SimplicialLDLT<SpMat> chol(A); // performs a Cholesky factorization of A

  Eigen::SparseVector<double> sp_b = b.sparseView(); // create sparse vector
  Eigen::VectorXd x = chol.solve(sp_b); // use the factorization to solve for the given right hand side

  for (int i = 1; i < nb_solve; i++) x = chol.solve(sp_b);

  return Rcpp::wrap(x);
}

/******************************************************************************/
