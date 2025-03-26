#include <RcppEigen.h>
using namespace Rcpp;

//' mat_mult
//'
//' Calculate matrix multiplication
//' @param A matrix
//' @param B matrix
//' @return A %*% B
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd mat_mult(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B) {
  return A * B ;  // Equivalent to A %*% B
}

