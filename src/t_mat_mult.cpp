#include <RcppEigen.h>
using namespace Rcpp;

// You can source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd t_mat_mult(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, const Eigen::MatrixXd& C) {
  return (A.transpose() * B) * C;  // Equivalent to (t(A) %*% B) %*% t(C)
}

