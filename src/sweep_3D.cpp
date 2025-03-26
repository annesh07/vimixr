#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Eigen;
using namespace Rcpp;

//' sweep_3D
//'
//' A C++ alternative of sweep() function from base R
//' @param A a 3D array
//' @param R a vector
//' @param dims dimensions in 3D
//' @param n_threads number of threads
//' @return sweep(A, 3, R, "*")
//' @export
// [[Rcpp::export]]
NumericVector sweep_3D(NumericVector A, NumericVector R, IntegerVector dims, int n_threads = 4) {
  int D = dims[0];    // First dimension
  int T0 = dims[2];   // Third dimension (number of slices)
  int D2 = D * D;     // Size of each slice

  // Map A and R to Eigen vectors
  Map<VectorXd> vec_A(Rcpp::as<Map<VectorXd>>(A));
  Map<VectorXd> vec_R(Rcpp::as<Map<VectorXd>>(R));

  // Clone A to avoid modifying input
  NumericVector result = clone(A);
  Map<VectorXd> vec_result(Rcpp::as<Map<VectorXd>>(result));

  // Parallelize over T0 slices
#pragma omp parallel for num_threads(n_threads)
  for (int t = 0; t < T0; t++) {
    int offset = t * D2;  // Start of slice t in flat memory
    vec_result.segment(offset, D2) *= vec_R[t];  // Scale the slice
  }

  // Preserve dimensions
  result.attr("dim") = dims;
  return result;
}
