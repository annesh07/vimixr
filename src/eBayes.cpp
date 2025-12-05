#include <RcppEigen.h>
using namespace Rcpp;

// --- global variables ---
Eigen::VectorXd Rk_g;
Eigen::MatrixXd PSk_g;
Eigen::MatrixXd Sk_g;

// --- initialization function ---
//' Initialize matrices for eBayes
//' 
//' @param Rk vector
//' @param PSk matrix
//' @param Sk matrix
//' @rdname eBayes
//' @export
// [[Rcpp::export]]
void eBayes_init(const Eigen::VectorXd& Rk,
                const Eigen::MatrixXd& PSk,
                const Eigen::MatrixXd& Sk) {
 Rk_g  = Rk;
 PSk_g = PSk;
 Sk_g  = Sk;
}
 
// --- run function ---
//' Compute eBayes for a_vec using global matrices
//'
//' @param a_vec vector
//' @return vector of results
//' @export
// [[Rcpp::export]]
Eigen::VectorXd eBayes_run(const Eigen::VectorXd& a_vec) {
  int n0 = a_vec.size();
  int T0 = Rk_g.size();
  int D  = PSk_g.cols();
  
  Eigen::VectorXd out(n0);
  
  // --- f0 computation for all i (vectorized over n0) ---
  // a_vec: n0 x 1, Rk_g: T0 x 1
  Eigen::MatrixXd f0_temp = a_vec.replicate(1, T0); // n0 x T0
  f0_temp = f0_temp.array() + Rk_g.transpose().replicate(n0, 1).array(); // ai + Rk_g
  Eigen::MatrixXd f0_num = (2 * f0_temp.array() + 1).array() * Rk_g.transpose().replicate(n0, 1).array();
  Eigen::MatrixXd f0_denom = f0_temp.array().square();
  Eigen::VectorXd out1 = 0.5 * D * f0_num.array().rowwise().sum() / f0_denom.array().rowwise().sum();
  
  // --- f1 computation (fully vectorized over n0 x T0 x D) ---
  // Repeat a_vec over T0 rows
  Eigen::MatrixXd a_mat = a_vec.replicate(1, T0);   // n0 x T0
  Eigen::MatrixXd Rk_mat = Rk_g.transpose().replicate(n0, 1); // n0 x T0
  Eigen::MatrixXd aiR = a_mat.array() * Rk_mat.array();       // n0 x T0
  
  // Expand aiR to n0 x T0 x D by broadcasting over columns
  Eigen::MatrixXd aiR_big = aiR.replicate(1, D);               // n0 x (T0*D)
  Eigen::MatrixXd PSk_big = PSk_g.transpose().replicate(n0, 1); // n0 x (T0*D)
  Eigen::MatrixXd Sk_big = Sk_g.transpose().replicate(n0, 1);   // n0 x (T0*D)
  
  Eigen::MatrixXd denom = (aiR_big.array() * 0 + 1.0 + Sk_big.array()).square(); // ai + Sk, but ai added in next step
  denom.array() += a_vec.replicate(1, T0*D).array(); // add ai
  denom = denom.array().square();                    // square
  
  Eigen::MatrixXd numer = aiR_big + PSk_big;
  
  // Sum each row
  Eigen::VectorXd out2 = (numer.array() / denom.array()).rowwise().sum();
  
  // Combine
  out = 0.5 * (out1 - out2);
  
  return out;
}