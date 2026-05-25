#include <Rcpp.h>
using namespace Rcpp;
//' sparse_cov_op
//'
//' Calculate the sum, squared sum and log sum of off-diagonal vector elements from the covariance array
//' @param x vector
//' @return a vector of sum, squared sum and log sum elements
//' @export
// [[Rcpp::export]]
 
 double sparse_cov_op(
     const Rcpp::NumericMatrix& X,      // n x p
     const Rcpp::NumericMatrix& P,      // n x K
     const Rcpp::NumericMatrix& inv_C0, // K x p
     const Rcpp::NumericMatrix& L1)     // K x p
 {
   const int n = X.nrow();
   const int p = X.ncol();
   const int K = P.ncol();
   double total = 0.0;
   
   for (int i = 0; i < n; i++) {
     for (int k = 0; k < K; k++) {
       const double pik = P(i, k);
       if (pik == 0.0) continue;
       double acc = 0.0;
       for (int j = 0; j < p; j++) {
         const double xij = X(i, j);
         acc += xij * (inv_C0(k,j) * (L1(k,j) - 0.5 * xij));
       }
       total += pik * acc;
     }
   }
   return total;
 }