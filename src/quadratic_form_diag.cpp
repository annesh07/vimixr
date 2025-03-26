#include <RcppEigen.h>
using namespace Rcpp;

//' quadratic_form_diag
//'
//' Calculate a combination of matrix multiplications
//' @param A matrix
//' @param B matrix
//' @return diag(A %*% B %*% t(A))
//' @export
// [[Rcpp::export]]
Eigen::VectorXd quadratic_form_diag(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B) {
  Eigen::MatrixXd AB = A * B;                          // Matrix product A %*% B
  Eigen::VectorXd diagResult = (AB.array() * A.array()).rowwise().sum(); // Element-wise product and row sums
  return diagResult;
}

