#include <Rcpp.h>
using namespace Rcpp;
//' lower_tri_stats
//'
//' Extract lower diagonal elements of a Matrix, and perform sum, squared sum and log sum
//' @param M matrix
//' @return a vector of sum, squared sum and log sum elements
//' @export
// [[Rcpp::export]]

Rcpp::NumericVector lower_tri_stats(const Rcpp::NumericMatrix& M) {
  int n = M.nrow();
  double s_x = 0.0, s_x2 = 0.0, s_logx = 0.0;
  
  for (int j = 0; j < n; j++) {
    const double* col = &M(0, j);   
    for (int i = j + 1; i < n; i++) {
      const double v = col[i];
      s_x    += v;
      s_x2   += v * v;
      s_logx -= std::log(v);
    }
  }
  return Rcpp::NumericVector::create(s_x, s_x2, s_logx);
}