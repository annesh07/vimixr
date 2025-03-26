#include <Rcpp.h>
using namespace Rcpp;
//' cum_clustprop_var Calculate the columnwise sum of rowwise cummulative probability for variance
//' @param P1 probability matrix
//' @return rowwise cummulative probability for variance
//' @export
// [[Rcpp::export]]

NumericVector cum_clustprop_var(const NumericMatrix P1) {
  int nrows = P1.nrow();
  int ncols = P1.ncol();
  NumericVector out(ncols, 0.0);

  // Process each row of P1
  for (int i = 0; i < nrows; i++) {
    // First, compute the total sum T for the row.
    double rowTotal = 0.0;
    for (int j = 0; j < ncols; j++) {
      rowTotal += P1(i, j);
    }
    // Then, compute the cumulative sum on the fly and update out[j].
    double cumulative = 0.0;
    for (int j = 0; j < ncols; j++) {
      cumulative += P1(i, j);
      // For column j, add contribution S(j)*(T - S(j))
      out[j] += cumulative * (rowTotal - cumulative);
    }
  }
  return out;
}
