#' Collapsed variational inference for non-parametric Bayesian mixture models
#'
#' @details The following models are supported in \code{vimixr}, listing their
#' required input arguments in \code{...} when calling \code{cvi_npmm()}: \itemize{
#'  \item \strong{Known covariance} \itemize{
#'        \item \emph{diagonal covariance} We need the followingadditional arguments: \describe{
#'            \item{\code{cov_data}: a scalar}{}
#'            }
#'          }
#'  \item \strong{Unknown covariance}
#' }
#'
#' @param X input data as a matrix
#' @param variational_params
#' @param prior_shape_alpha
#' @param prior_rate_alpha
#' @param post_shape_alpha
#' @param post_rate_alpha
#' @param prior_mean_eta
#' @param post_mean_eta
#' @param log_prob_matrix
#' @param maxit maximum number of iterations. Default is 100
#' @param fixed_variance covariance matrix of the data is considered known (fixed)
#' or unknown. Default is FALSE
#' @param covariance_type covariance matrix is considered diagonal or full.
#' Default is 'full'
#' @param cluster_specific_covariance covariance matrix is specific to a cluster
#' allocation or it is same over all cluster choices. Default is TRUE
#' @param variance_prior_type For unknown and full covariance matrix, choice of
#' matrix prior is either Inverse-Wishart ('IW') or Cholesky-decomposed
#' ('decomposed'). For unknown, full and cluster-specific covariance matrix,
#' choice of matrix prior is either Inverse-Wishart ('IW'), element-wise Gamma
#' and Laplace distributed ('sparse') or element-wise Gamma and Normal
#' distributed ('off-diagonal normal')
#' @param ... additional paremeters for specific models. See Details below.
#'
#' @returns Posterior DP concentration parameter ('alpha'),
#' Posterior number of clusters ('Cluster number'),
#' Posterior proportions of clusters ('Cluster Proportion'),
#' Posterior logarithm of cluster allocation matrix ('log Probability matrix')
#' and Optimisation of the ELBO function ('ELBO')
#'
#' @importFrom Rfast rowsums colsums spdinv Crossprod Tcrossprod mat.mult
#' Diag.fill Diag.matrix
#'
#' @export
#'
#' @examples
#'
#' X <- rbind(matrix(rnorm(100, m=0, sd=0.5), ncol=2),
#'            matrix(rnorm(100, m=3, sd=0.5), ncol=2)
#'            )
#' cvi_npmm(X, fixed_variance=TRUE)
#'
#'
cvi_npmm <- function(X, variational_params,
                     prior_shape_alpha, prior_rate_alpha,
                     post_shape_alpha, post_rate_alpha,
                     prior_mean_eta, post_mean_eta,
                     log_prob_matrix,
                     maxit = 100,
                     fixed_variance=FALSE, covariance_type="full",
                     cluster_specific_covariance=TRUE,
                     variance_prior_type=c("IW", "decomposed", "sparse",
                                           "off-diagonal normal"),
                     ...){

  N <- nrow(X) #number of samples
  D <- ncol(X) # number of variables
  T0 <- variational_params

  params <- list()

  params$prior_mean_eta <- prior_mean_eta
  params$prior_shape_alpha <- prior_shape_alpha
  params$prior_rate_alpha <- prior_rate_alpha
  params$post_shape_alpha <- post_shape_alpha
  params$post_rate_alpha <- post_rate_alpha
  params$post_mean_eta <- post_mean_eta

  params$log_prob_matrix <- log_prob_matrix
  params$P <- exp(log_prob_matrix)
  RP <- colsums(params$P)

  #updating the parameter list based on the conditions
  if(covariance_type == "diagonal") {

    if(fixed_variance) {
      params$prior_precision_scalar_eta <- prior_precision_scalar_eta
      params$post_precision_scalar_eta <- post_precision_scalar_eta
      params$cov_data <- cov_data
      params_check(params, fixed_variance, covariance_type,
                   cluster_specific_covariance,
                   variance_prior_type)

      C00 <- diag(D)/prior_precision_scalar_eta #covariance of DP mean parameters
      inverts[["inv_C0"]] <- Rfast::spdinv(cov_data)
      inverts[["inv_C00"]] <- Rfast::spdinv(C00)

    } else {
      params$prior_shape_scalar_cov <- prior_shape_scalar_cov
      params$prior_rate_scalar_cov <- prior_rate_scalar_cov
      params$post_shape_scalar_cov <- post_shape_scalar_cov
      params$post_rate_scalar_cov <- post_rate_scalar_cov
      params$post_precision_scalar_eta <- post_precision_scalar_eta
      params$prior_precision_scalar_eta <- prior_precision_scalar_eta
      params_check(params, fixed_variance, covariance_type,
                   cluster_specific_covariance,
                   variance_prior_type)

      C00 <- diag(D)/prior_precision_scalar_eta #covariance of DP mean parameters
      inverts[["inv_C00"]] <- Rfast::spdinv(C00)

    }

  } else if(covariance_type == "full") {

    if(fixed_variance) {
      params$post_cov_eta <- post_cov_eta
      params$cov_data <- cov_data
      params$prior_cov_eta <- prior_cov_eta
      params_check(params, fixed_variance, covariance_type,
                   cluster_specific_covariance,
                   variance_prior_type)

      inverts[["inv_C0"]] <- Rfast::spdinv(cov_data)
      inverts[["inv_C00"]] <- Rfast::spdinv(prior_cov_eta)


    } else {
      if(!cluster_specific_covariance) {
        if(variance_prior_type == "IW"){
          params$prior_df_cov <- prior_df_cov
          params$prior_scale_cov <- prior_scale_cov
          params$post_df_cov <- post_df_cov
          params$post_scale_cov <- post_scale_cov
          params$post_cov_eta <- post_cov_eta
          params$prior_cov_eta <- prior_cov_eta
          params_check(params, fixed_variance, covariance_type,
                       cluster_specific_covariance,
                       variance_prior_type)

          inverts[["inv_V0"]] <- Rfast::spdinv(prior_scale_cov)
          inverts[["inv_C00"]] <- Rfast::spdinv(prior_cov_eta)


        } else if (variance_prior_type == "decomposed"){
          params$prior_scale_diag_decomp <- prior_scale_diag_decomp
          params$prior_rate_diag_decomp <- prior_rate_diag_decomp
          params$prior_mean_offdiag_decomp <- prior_mean_offdiag_decomp
          params$prior_var_offdiag_decomp <- prior_var_offdiag_decomp
          params$post_scale_diag_decomp <- post_scale_diag_decomp
          params$post_rate_diag_decomp <- post_rate_diag_decomp
          params$post_mean_offdiag_decomp <- post_mean_offdiag_decomp
          params$post_var_offdiag_decomp <- post_var_offdiag_decomp
          params$post_cov_eta <- post_cov_eta
          params$prior_cov_eta <- prior_cov_eta
          params_check(params, fixed_variance, covariance_type,
                       cluster_specific_covariance,
                       variance_prior_type)

          inverts[["inv_C00"]] <- Rfast::spdinv(prior_cov_eta)


        } else {
          stop("'variance_prior_type' can only be either 'IW' or 'decomposed'
               when 'cluster_specific_covariance' is FALSE")
        }

      }else{
        if(variance_prior_type == "IW"){
          params$prior_df_cs_cov <- prior_df_cs_cov
          params$prior_scale_cs_cov <- prior_scale_cs_cov
          params$post_df_cs_cov <- post_df_cs_cov
          params$post_scale_cs_cov <- post_scale_cs_cov
          params$scaling_cov_eta <- scaling_cov_eta
          params_check(params, fixed_variance, covariance_type,
                       cluster_specific_covariance,
                       variance_prior_type)


        } else if (variance_prior_type == "sparse"){
          params$prior_scale_d_cs_cov <- prior_scale_d_cs_cov
          params$prior_rate_d_cs_cov <- prior_rate_d_cs_cov
          params$prior_var_offd_cs_cov <- prior_var_offd_cs_cov
          params$post_scale_d_cs_cov <- post_scale_d_cs_cov
          params$post_rate_d_cs_cov <- post_rate_d_cs_cov
          params$post_var_offd_cs_cov <- post_var_offd_cs_cov
          params$scaling_cov_eta <- scaling_cov_eta
          params_check(params, fixed_variance, covariance_type,
                       cluster_specific_covariance,
                       variance_prior_type)


        } else if (variance_prior_type == "off-diagonal normal"){
          params$prior_scale_d_cs_cov <- prior_scale_d_cs_cov
          params$prior_rate_d_cs_cov <- prior_rate_d_cs_cov
          params$post_scale_d_cs_cov <- post_scale_d_cs_cov
          params$post_rate_d_cs_cov <- post_rate_d_cs_cov
          params$post_mean_offd_cs_cov <- post_mean_offd_cs_cov
          params$scaling_cov_eta <- scaling_cov_eta
          params_check(params, fixed_variance, covariance_type,
                       cluster_specific_covariance,
                       variance_prior_type)


        } else {
          stop("'variance_prior_type' can only be either 'IW' or 'decomposed'
               when 'cluster_specific_covariance' is TRUE")
        }
      }

    }
  } else {
    stop("covariance_type can only be either 'diagonal' or 'full'.")
  }

  #store the output of ELBO function for every iteration of updates
  elbo_values <- list()
  elbo_values[[1]] <- ELBO_function(fixed_variance, covariance_type,
                                    cluster_specific_covariance,
                                    variance_prior_type, X, inverts, params)

  for (m in 1:maxit){

    updated_params <- CVI_update_function(fixed_variance, covariance_type,
                                          cluster_specific_covariance,
                                          variance_prior_type,
                                          X, inverts, params)

    params <- updated_params

    elbo_values[[m+1]] <- ELBO_function(fixed_variance, covariance_type,
                                        cluster_specific_covariance,
                                        variance_prior_type, X, inverts, params)

    if (abs(sum(elbo_values[[m]]) - sum(elbo_values[[m + 1]])) < 0.000001 ){
      break
    }
    message("outer loop: ", m,"\n", elbo_values[[m + 1]], '\n', sep="")
  }
  W1 <- params$post_shape_alpha
  W2 <- params$post_rate_alpha
  Plog <- params$log_prob_matrix

  alpha0 <- W1/W2 #posterior concentration parameter
  clustering <- apply(Plog, MARGIN = 1, FUN=which.max)
  clust <- table(clustering) #clusters with proportions
  clustnum <- length(unique(clustering)) #number of clusters
  posterior <- list("alpha" = alpha0, "Cluster number" = clustnum,
                    "Cluster Proportion" = clust,
                    "log Probability matrix" = Plog)
  optimisation <- list("ELBO" = elbo_values)

  output <-  list("posterior" = posterior, "optimisation" = optimisation)
  class(output) <- "CVIoutput"

  return(output)
}

