#' Generate random log Probability matrix if not provided 
#'
#' @param N rows of the data matrix
#' @param T0 variational clusters
#' @param seed0 seed for generating log Probability matrix
#' 
#' @return No return value, called for side effects.
#' 
#' @importFrom stats rgamma
#' 
#' @export
generate_log_prob <- function(N, T0, seed0) {
  set.seed(seed0)
  probs <- matrix(stats::rgamma(N * T0, shape = 1), nrow = N, ncol = T0)
  probs <- probs / rowSums(probs)
  log(probs)
}