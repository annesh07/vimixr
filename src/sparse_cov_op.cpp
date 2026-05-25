#include <Rcpp.h>
using namespace Rcpp;
//' sparse_cov_op
//'
//' Calculate the sum, squared sum and log sum of off-diagonal vector elements from the covariance array
//' @param X data matrix
//' @param P probability allocation matrix
//' @param inv_C0 matrix corresponding diagonal elements of the cluster precision matrices
//' @param L1 cluster mean matrix
//' @return likelihood term calculation in elbo
//' @export
// [[Rcpp::export]]
 
 double sparse_cov_op(
     const Rcpp::NumericMatrix& X,      
     const Rcpp::NumericMatrix& P,      
     const Rcpp::NumericMatrix& inv_C0, 
     const Rcpp::NumericMatrix& L1)     
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