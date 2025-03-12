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
Eigen::MatrixXd mat_mult_t(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, const Eigen::MatrixXd& C) {
  return A * (B * C.transpose());  // Equivalent to A %*% (B %*% t(C))
}

