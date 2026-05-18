#include <RcppEigen.h>
using namespace Rcpp;

//' mat_mult
//'
//' Calculate matrix multiplication with optional transposition.
//' @param A matrix or vector
//' @param B matrix or vector
//' @param transpose_A transpose A before multiplying
//' @param transpose_B transpose B before multiplying
//' @return A %*% B (or variant), as vector if either input was a vector
//' @export
// [[Rcpp::export]]
SEXP mat_mult(SEXP A, SEXP B,
             bool transpose_A = false,
             bool transpose_B = false) {
 
 bool a_is_vec = !Rf_isMatrix(A);
 bool b_is_vec = !Rf_isMatrix(B);
 
 auto toMatrix = [](SEXP x, bool is_vec) -> Eigen::MatrixXd {
   if (!is_vec) return Rcpp::as<Eigen::MatrixXd>(x);
   Eigen::VectorXd v = Rcpp::as<Eigen::VectorXd>(x);
   return v; 
 };
 
 Eigen::MatrixXd matA = toMatrix(A, a_is_vec);
 Eigen::MatrixXd matB = toMatrix(B, b_is_vec);
 
 Eigen::MatrixXd result(
     transpose_A ? matA.cols() : matA.rows(),
     transpose_B ? matB.rows() : matB.cols()
 );
 
 if      ( transpose_A && !transpose_B) result.noalias() = matA.transpose() * matB;
 else if (!transpose_A &&  transpose_B) result.noalias() = matA * matB.transpose();
 else if ( transpose_A &&  transpose_B) result.noalias() = matA.transpose() * matB.transpose();
 else                                   result.noalias() = matA * matB;
 
 if (a_is_vec || b_is_vec) {
   return Rcpp::wrap(Eigen::VectorXd(result.reshaped()));
 }
 return Rcpp::wrap(result);
}