params_check <- function(params, fixed_variance=FALSE, covariance_type="diagonal",
                         cluster_specific_covariance=TRUE,
                         variance_prior_type=c("IW", "decomposed", "sparse",
                                               "off-diagonal normal")){

  if(covariance_type == "diagonal") {

    if(fixed_variance) {
      if(is.null(params$prior_precision_scalar_eta) &&
         is.null(params$post_precision_scalar_eta) &&
         is.null(params$cov_data)
      ){
        stop("fixed_variance=TRUE and covariance_type == 'diagonal' require the following arguments to be non-null:
             'prior_precision_scalar_eta', 'post_precision_scalar_eta', 'cov_data'")
      }

    } else {
      params$prior_shape_scalar_cov <- prior_shape_scalar_cov
      params$prior_rate_scalar_cov <- prior_rate_scalar_cov
      params$post_shape_scalar_cov <- post_shape_scalar_cov
      params$post_rate_scalar_cov <- post_rate_scalar_cov
      params$post_precision_scalar_eta <- post_precision_scalar_eta
      params$prior_precision_scalar_eta <- prior_precision_scalar_eta

      C00 <- diag(D)/prior_precision_scalar_eta #covariance of DP mean parameters
      inverts[["inv_C00"]] <- Rfast::spdinv(C00)

    }

  } else if(covariance_type == "full") {

    if(fixed_variance) {
      params$post_cov_eta <- post_cov_eta
      params$cov_data <- cov_data
      params$prior_cov_eta <- prior_cov_eta

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


        } else if (variance_prior_type == "sparse"){
          params$prior_scale_d_cs_cov <- prior_scale_d_cs_cov
          params$prior_rate_d_cs_cov <- prior_rate_d_cs_cov
          params$prior_var_offd_cs_cov <- prior_var_offd_cs_cov
          params$post_scale_d_cs_cov <- post_scale_d_cs_cov
          params$post_rate_d_cs_cov <- post_rate_d_cs_cov
          params$post_var_offd_cs_cov <- post_var_offd_cs_cov
          params$scaling_cov_eta <- scaling_cov_eta


        } else if (variance_prior_type == "off-diagonal normal"){
          params$prior_scale_d_cs_cov <- prior_scale_d_cs_cov
          params$prior_rate_d_cs_cov <- prior_rate_d_cs_cov
          params$post_scale_d_cs_cov <- post_scale_d_cs_cov
          params$post_rate_d_cs_cov <- post_rate_d_cs_cov
          params$post_mean_offd_cs_cov <- post_mean_offd_cs_cov
          params$scaling_cov_eta <- scaling_cov_eta


        } else {
          stop("'variance_prior_type' can only be either 'IW' or 'decomposed'
               when 'cluster_specific_covariance' is TRUE")
        }
      }

    }
  } else {
    stop("covariance_type can only be either 'diagonal' or 'full'.")
  }



}
