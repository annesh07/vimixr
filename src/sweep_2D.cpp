#include <RcppEigen.h>
#include <string>
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
Eigen::MatrixXd sweep_2D(Eigen::MatrixXd& A, const Eigen::VectorXd& B,
                         const std::string& operation, const int margin) {

  if (margin == 1) {  // Column-wise operation
    if (operation == "+") {
      A.array().colwise() += B.array();  // Element-wise addition
    } else if (operation == "-") {
      A.array().colwise() -= B.array();  // Element-wise subtraction
    } else if (operation == "*") {
      A.array().colwise() *= B.array();  // Element-wise multiplication
    } else if (operation == "/") {
      A.array().colwise() /= B.array();  // Element-wise division
    } else {
      stop("Invalid operation. Use '+', '-', '*', or '/'");
    }
  }
  else if (margin == 2) {  // Row-wise operation
    if (operation == "+") {
      A.array().rowwise() += B.transpose().array();  // Element-wise addition
    } else if (operation == "-") {
      A.array().rowwise() -= B.transpose().array();  // Element-wise subtraction
    } else if (operation == "*") {
      A.array().rowwise() *= B.transpose().array();  // Element-wise multiplication
    } else if (operation == "/") {
      A.array().rowwise() /= B.transpose().array();  // Element-wise division
    } else {
      stop("Invalid operation. Use '+', '-', '*', or '/'");
    }
  }
  else {
    stop("Invalid margin. Use 1 for column-wise or 2 for row-wise.");
  }

  return A;  // Return the modified matrix A
}
