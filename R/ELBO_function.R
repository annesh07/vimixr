#' General ELBO function
#'
#' @param fixed_variance whether the covariance is fixed or estimated.
#' Default is \code{FALSE} which means it is estimated.
#' @param covariance_type The assumed type of the covariance matrix.
#' Can be either \code{"diagonal"} if it is the identify multiplied by a scalar,
#' or \code{"full"} for a fully unspecified covariance matrix.
#' @param cluster_specific_covariance whether the the covariance is shared across
#' estimated clusters or is cluster specific. Default is \code{TRUE} which means it is cluster specific.
#' @param variance_prior_type character string specifying the type of prior distribution
#' for the covariance when cluster_specific_covariance is \code{TRUE}.
#' Can be either \code{"IW"} or \code{"decomposed"} if \code{cluster_specific_covariance} is \code{FALSE},
#' and can be either \code{"IW"}, \code{"sparse"} or \code{"off-diagonal normal"} otherwise.
#' @param X the data matrix
#' @param inverts a list of inverses
#' @param params a list of required arguments
#'
#' @returns ELBO values
#'
#' @export
#'
#' @examples
ELBO_function <- function(fixed_variance = FALSE, covariance_type = "diagonal",
                          cluster_specific_covariance = TRUE,
                          variance_prior_type = c("IW", "decomposed", "sparse",
                                                  "off-diagonal normal"),
                          X,
                          inverts,
                          params){
  N <- params$N
  D <- params$D
  T0 <- params$T0
  L1 <- params$post_mean_eta
  Mu0 <- params$prior_mean_eta
  s1 <- params$prior_shape_alpha #shape1 parameter for alpha prior
  s2 <- params$prior_rate_alpha #shape2 parameter for alpha prior
  W1 <- params$post_shape_alpha #posterior shape1 parameter for alpha
  W2 <- params$post_rate_alpha #posterior shape2 parameter for alpha

  Plog <- params[["log_prob_matrix"]] #log of posterior probability allocation matrix
  P <- params[["P"]]
  RP <- Rfast::colsums(P)

  #expectation of alpha prior
  e_alpha <- s1*log(s2) - lgamma(s1) + (s1 - 1)*(digamma(W1) - log(W2)) -
    s2 * (W1 / W2)

  #expectation of latent probability allocation prior
  cp <- RP #cluster proportions
  v_cp <- Rfast::colsums(P*(1 - P)) #variance of the cluster proportions
  ccp <- cum_clustprop(P) #cumulative cluster proportions
  v_ccp <- cum_clustprop_var(P) #variance of the cumulative cluster proportions
  e_indiv_alloc <- lgamma(1 + cp) + 0.5 * trigamma(1 + cp) * v_cp +
    lgamma(W1 / W2 + ccp) + 0.5 * trigamma(W1 / W2 + ccp) * ((W1 / W2^2) + v_ccp) -
    lgamma(1 + W1 / W2 + cp + ccp) -
    0.5*trigamma(1 + W1 / W2 + cp + ccp) * (W1 / W2^2 + v_cp + v_ccp)
  e_alloc <- T0 * (digamma(W1) - log(W2)) + sum(e_indiv_alloc)

  #Variational expectation of alpha & latent allocations
  e_alpha_post <- W1 * log(W2) - lgamma(W1) + (W1 - 1)*(-log(W2) + digamma(W1)) - W1
  e_alloc_post <- sum(exp(Plog) * Plog)

  if(covariance_type == "diagonal") {

    if(fixed_variance) {

      fixed_diagonal_elbo <- elbo_fixed_diagonal(X, inverts, params)
      fixed_diagonal_elbo[["me_var"]] <- fixed_diagonal_elbo[["me_var"]] -
        e_alpha_post - e_alloc_post
      out <- c("e_alpha" = e_alpha,  "e_alloc" = e_alloc,
               fixed_diagonal_elbo)

    } else {

      varied_diagonal_elbo <- elbo_varied_diagonal(X, inverts, params)
      varied_diagonal_elbo[["me_var"]] <- varied_diagonal_elbo[["me_var"]] -
        e_alpha_post - e_alloc_post
      out <- c("e_alpha" = e_alpha,  "e_alloc" = e_alloc,
               varied_diagonal_elbo)

    }

  } else if(covariance_type == "full") {

    if(fixed_variance) {

      fixed_full_elbo <- elbo_fixed_full(X, inverts, params)
      fixed_full_elbo[["me_var"]] <- fixed_full_elbo[["me_var"]] -
        e_alpha_post - e_alloc_post
      out <- c("e_alpha" = e_alpha,  "e_alloc" = e_alloc,
               fixed_full_elbo)

    } else {
      if(!cluster_specific_covariance) {
        if(variance_prior_type == "IW"){

          varied_IW_full_elbo <- elbo_varied_IW_full(X, inverts, params)
          varied_IW_full_elbo[["me_var"]] <- varied_IW_full_elbo[["me_var"]] -
            e_alpha_post - e_alloc_post
          out <- c("e_alpha" = e_alpha,  "e_alloc" = e_alloc,
                   varied_IW_full_elbo)

        } else if (variance_prior_type == "decomposed"){

          varied_decomposed_full_elbo <- elbo_varied_decomposed_full(X, inverts,
                                                                     params)
          varied_decomposed_full_elbo[["me_var"]] <- varied_decomposed_full_elbo[["me_var"]] -
            e_alpha_post - e_alloc_post
          out <- c("e_alpha" = e_alpha,  "e_alloc" = e_alloc,
                   varied_decomposed_full_elbo)

        } else {
          stop("'variance_prior_type' can only be either 'IW' or 'decomposed'
               when 'cluster_specific_covariance' is FALSE")
        }

      }else{
        if(variance_prior_type == "IW"){

          cs_IW_elbo <- elbo_cs_IW(X, inverts, params)
          cs_IW_elbo[["me_var"]] <- cs_IW_elbo[["me_var"]] -
            e_alpha_post - e_alloc_post
          out <- c("e_alpha" = e_alpha,  "e_alloc" = e_alloc,
                   cs_IW_elbo)

        } else if (variance_prior_type == "sparse"){

          cs_sparse_elbo <- elbo_cs_sparse(X, inverts, params)
          cs_sparse_elbo[["me_var"]] <- cs_sparse_elbo[["me_var"]] -
            e_alpha_post - e_alloc_post
          out <- c("e_alpha" = e_alpha,  "e_alloc" = e_alloc,
                   cs_sparse_elbo)

        } else if (variance_prior_type == "off-diagonal normal"){

          cs_offd_normal_elbo <- elbo_cs_offd_normal(X, inverts, params)
          cs_offd_normal_elbo[["me_var"]] <- cs_offd_normal_elbo[["me_var"]] -
            e_alpha_post - e_alloc_post
          out <- c("e_alpha" = e_alpha,  "e_alloc" = e_alloc,
                   cs_offd_normal_elbo)

        } else {
          stop("'variance_prior_type' can only be either 'IW' or 'decomposed'
               when 'cluster_specific_covariance' is TRUE")
        }
      }

    }
  } else {
    stop("covariance_type can only be either 'diagonal' or 'full'.")
  }

  return("ELBO" = out)
}
