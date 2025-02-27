#' Title
#'
#' @param Data
#' @param variational_params
#' @param prior_shape_alpha
#' @param prior_rate_alpha
#' @param post_shape_alpha
#' @param post_rate_alpha
#' @param prior_mean_eta
#' @param post_mean_eta
#' @param log_prob_matrix
#' @param maxit
#' @param type
#' @param params
#'
#' @returns
#' @export
#'
#' @examples
CollapsedVI <- function(Data, variational_params,
                        prior_shape_alpha, prior_rate_alpha,
                        post_shape_alpha, post_rate_alpha,
                        prior_mean_eta, post_mean_eta,
                        log_prob_matrix,
                        prior_scale_d_cs_cov=NULL,
                        maxit,
                        fixed_variance=FALSE, covariance_type="diagonal",
                        cluster_specific_covariance=TRUE,
                        variance_prior_type=c("IW", "decomposed", "sparse", "off-diagonal normal")){
  library("Rfast")
  X <- Data
  N <- dim(X)[1]
  D <- dim(X)[2]
  T0 <- variational_params

  params <- list()

  params$prior_mean_eta <- prior_mean_eta
  params$prior_shape_alpha <- prior_shape_alpha
  params$prior_rate_alpha <- prior_rate_alpha
  params$post_shape_alpha <- post_shape_alpha
  params$post_rate_alpha <- post_rate_alpha
  params$post_mean_eta <- post_mean_eta

  params$Plog <- log_prob_matrix
  params$P <- exp(Plog)
  RP <- colsums(P)

  C0 <- params$cov_data
  L20 <- params$prior_precision_scalar_eta
  L2 <- params$post_precision_scalar_eta

  if(fixed_variance){

  }else{
    if(cluster_specific_covariance){
      if(variance_prior_type == "off-diagonal normal"){
        stopifnot(!is.null(prior_scale_d_cs_cov))
        params$prior_scale_d_cs_cov <- prior_scale_d_cs_cov
      }
    }
  }
  # computing inverts for the different cases
  C00 <- diag(D)/L20       #C00; covariance matrix for eta_i's
  inv_C00 <- spdinv(C00)    #inverse of C00
  inv_C0 <- spdinv(C0)      #inverse of C0

  #store the output of ELBO function for every iteration of updates
  elbo_values <- list()
  elbo_values[[1]] <- ELBO_function(fixed_variance, covariance_type,
                                    cluster_specific_covariance,
                                    variance_prior_type, X, inverts, params)

  for (m in 1:maxit){

    updated_params <- params_update(fixed_variance, covariance_type,
                                       cluster_specific_covariance,
                                       variance_prior_type, params)

    params <- updated_params

    elbo_values[m+1] <- ELBO_function(fixed_variance, covariance_type,
                                     cluster_specific_covariance,
                                     variance_prior_type, X, inverts, params)
  }



}

}
