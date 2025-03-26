#include <RcppEigen.h>
using namespace Rcpp;

//' t_mat_mult
//'
//' Calculate a combination of matrix multiplications
//' @param A matrix
//' @param B matrix
//' @param C matrix
//' @return t(A) %*% B %*% C
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd t_mat_mult(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, const Eigen::MatrixXd& C) {
  return (A.transpose() * B) * C;  // Equivalent to (t(A) %*% B) %*% t(C)
}

