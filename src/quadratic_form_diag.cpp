#include <RcppEigen.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
//' @export
// [[Rcpp::export]]
Eigen::VectorXd quadratic_form_diag(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B) {
  Eigen::MatrixXd AB = A * B;                          // Matrix product A %*% B
  Eigen::VectorXd diagResult = (AB.array() * A.array()).rowwise().sum(); // Element-wise product and row sums
  return diagResult;
}

