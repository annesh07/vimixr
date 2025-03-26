#include <RcppEigen.h>
using namespace Rcpp;

//' mat_mult_t
//'
//' Calculate a combination of matrix multiplications
//' @param A matrix
//' @param B matrix
//' @param C matrix
//' @return A %*% B %*% t(C)
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd mat_mult_t(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, const Eigen::MatrixXd& C) {
  return A * (B * C.transpose());  // Equivalent to A %*% (B %*% t(C))
}

