#' Generate random log Probability matrix if not provided 
#'
#' @param N rows of the data matrix
#' @param T0 variational clusters
#' @param seed seed for generating log Probability matrix
#' 
#' @importFrom stats rgamma
#' 
#' @export
generate_log_prob <- function(N, T0, seed) {
  set.seed(seed)
  probs <- matrix(stats::rgamma(N * T0, shape = 1), nrow = N, ncol = T0)
  probs <- probs / rowSums(probs)
  log(probs)
}

#' Generate random Sparse DPMM hyperparameters for cluster-specific covariance 
#' matrices
#'
#' @param T0 variational clusters
#' @param D dimensions of the data matrix
#' @param seed seed for generating log Probability matrix
#'
#' @importFrom stats runif
#' 
#' @export
generate_cs_priors <- function(T0, D, seed) {
  set.seed(seed)
  prior_shape_d_cs_cov <- matrix(stats::runif(T0, 0.001, D), nrow = 1)  
  prior_rate_d_cs_cov <- matrix(stats::runif(T0 * D, 0.001, D), nrow = T0, ncol = D)  
  list(prior_shape_d_cs_cov = prior_shape_d_cs_cov,
       prior_rate_d_cs_cov = prior_rate_d_cs_cov)
}