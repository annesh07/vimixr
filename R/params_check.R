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
        stop("fixed_variance = TRUE and covariance_type = 'diagonal' require the
              following arguments to be non-null: 'prior_precision_scalar_eta',
             'post_precision_scalar_eta', 'cov_data'")
      }

    } else {
      if(is.null(params$prior_shape_scalar_cov) &&
         is.null(params$prior_rate_scalar_cov) &&
         is.null(params$post_shape_scalar_cov) &&
         is.null(params$post_rate_scalar_cov) &&
         is.null(params$post_precision_scalar_eta) &&
         is.null(params$prior_precision_scalar_eta)
      ){
        stop("fixed_variance = FALSE and covariance_type = 'diagonal' require the
              following arguments to be non-null: 'prior_shape_scalar_cov',
             'prior_rate_scalar_cov', 'post_shape_scalar_cov',
             'post_rate_scalar_cov', 'post_precision_scalar_eta',
             'prior_precision_scalar_eta'")
      }
    }

  } else if(covariance_type == "full") {

    if(fixed_variance) {
      if(is.null(params$post_cov_eta) &&
         is.null(params$cov_data) &&
         is.null(params$prior_cov_eta)
      ){
        stop("fixed_variance = TRUE and covariance_type = 'full' require the
             following argument to be non-null: 'post_cov_eta', 'cov_data',
             'prior_cov_eta'")
      }

    } else {
      if(!cluster_specific_covariance) {
        if(variance_prior_type == "IW"){
          if(is.null(params$prior_df_cov) &&
             is.null(params$prior_scale_cov) &&
             is.null(params$post_df_cov) &&
             is.null(params$post_scale_cov) &&
             is.null(params$post_cov_eta) &&
             is.null(params$prior_cov_eta)
          ){
            stop("fixed_variance = FALSE, covariance_type = 'full',
                 cluster_specific_covariance = FALSE and
                 variance_prior_type = 'IW' require the following argument
                 to be non-null: 'prior_df_cov', 'prior_scale_cov',
                 'post_df_cov', 'post_scale_cov', 'post_cov_eta',
                 'prior_cov_eta'")
          }

        } else if (variance_prior_type == "decomposed"){
          if(is.null(params$prior_shape_diag_decomp) &&
             is.null(params$prior_rate_diag_decomp) &&
             is.null(params$prior_mean_offdiag_decomp) &&
             is.null(params$prior_var_offdiag_decomp) &&
             is.null(params$post_shape_diag_decomp) &&
             is.null(params$post_rate_diag_decomp) &&
             is.null(params$post_mean_offdiag_decomp) &&
             is.null(params$post_var_offdiag_decomp)  &&
             is.null(params$post_cov_eta) &&
             is.null(params$prior_cov_eta)
          ){
            stop("fixed_variance = FALSE, covariance_type = 'full',
                 cluster_specific_covariance = FALSE and
                 variance_prior_type = 'decomposed' require the following
                 argument to be non-null: 'prior_shape_diag_decomp',
                 'prior_rate_diag_decomp', 'prior_mean_offdiag_decomp',
                 'prior_var_offdiag_decomp', 'post_shape_diag_decomp',
                 'post_rate_diag_decomp', 'post_mean_offdiag_decomp',
                 'post_var_offdiag_decomp', 'post_cov_eta', 'prior_cov_eta'")
          }


        } else {
          stop("'variance_prior_type' can only be either 'IW' or 'decomposed'
               when 'cluster_specific_covariance' is FALSE")
        }

      }else{
        if(variance_prior_type == "IW"){
          if(is.null(params$prior_df_cs_cov) &&
             is.null(params$prior_scale_cs_cov) &&
             is.null(params$post_df_cs_cov) &&
             is.null(params$post_scale_cs_cov) &&
             is.null(params$scaling_cov_eta)
          ){
            stop("fixed_variance = FALSE, covariance_type = 'full',
                 cluster_specific_covariance = TRUE and
                 variance_prior_type = 'IW' require the following
                 argument to be non-null: 'prior_df_cs_cov',
                 'prior_scale_cs_cov', 'post_df_cs_cov',
                 'post_scale_cs_cov', 'scaling_cov_eta'")
          }

        } else if (variance_prior_type == "sparse"){
          if(is.null(params$prior_shape_d_cs_cov) &&
             is.null(params$prior_rate_d_cs_cov) &&
             is.null(params$prior_var_offd_cs_cov) &&
             is.null(params$post_shape_d_cs_cov) &&
             is.null(params$post_rate_d_cs_cov) &&
             is.null(params$post_var_offd_cs_cov) &&
             is.null(params$scaling_cov_eta)
          ){
            stop("fixed_variance = FALSE, covariance_type = 'full',
                 cluster_specific_covariance = TRUE and
                 variance_prior_type = 'sparse' require the following
                 argument to be non-null: 'prior_shape_d_cs_cov',
                 'prior_rate_d_cs_cov', 'prior_var_offd_cs_cov',
                 'post_shape_d_cs_cov', 'post_rate_d_cs_cov',
                 'post_var_offd_cs_cov', 'scaling_cov_eta'")
          }

        } else if (variance_prior_type == "off-diagonal normal"){
          if(is.null(params$prior_shape_d_cs_cov) &&
             is.null(params$prior_rate_d_cs_cov) &&
             is.null(params$post_shape_d_cs_cov) &&
             is.null(params$post_rate_d_cs_cov) &&
             is.null(params$post_mean_offd_cs_cov) &&
             is.null(params$scaling_cov_eta)
          ){
            stop("fixed_variance = FALSE, covariance_type = 'full',
                 cluster_specific_covariance = TRUE and
                 variance_prior_type = 'off-diagonal normal'
                 require the following argument to be non-null:
                 'prior_shape_d_cs_cov', 'prior_rate_d_cs_cov',
                 'post_shape_d_cs_cov', 'post_rate_d_cs_cov',
                 'post_mean_offd_cs_cov', 'scaling_cov_eta'")
          }

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
