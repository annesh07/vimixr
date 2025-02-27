CollapsedVI <- function(Data, variational_params,
                        prior_shape_alpha, prior_rate_alpha,
                        post_shape_alpha, post_rate_alpha,
                        prior_mean_eta, post_mean_eta,
                        log_prob_matrix, maxit, type, params){
  library("Rfast")
  X <- Data
  N <- dim(X)[1]
  D <- dim(X)[2]
  T0 <- variational_params
  s1 <- prior_shape_alpha
  s2 <- prior_rate_alpha
  Mu0 <- prior_mean_eta

  W1 <- post_shape_alpha
  W2 <- post_rate_alpha
  L1 <- post_mean_eta

  Plog <- log_prob_matrix
  P <- exp(Plog)
  RP <- colsums(P)

  if (type == "fixed-scalar"){
    C0 <- params$cov_data
    L20 <- params$prior_precision_scalar_eta
    L2 <- params$post_precision_scalar_eta

    C00 <- diag(D)/L20       #C00; covariance matrix for eta_i's
    inv_C00 <- spdinv(C00)    #inverse of C00
    inv_C0 <- spdinv(C0)      #inverse of C0

    #store the output of ELBO function for every iteration of updates
    f <- list()
    f[[1]] <- ELBO_function(type, params)

    for (m in 1:maxit){

    }



  }

}
