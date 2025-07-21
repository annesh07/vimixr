#' Collapsed variational inference for non-parametric Bayesian mixture models
#'
#' @details The following models are supported in \code{vimixr}, listing their
#' required input arguments in \code{...} when calling \code{cvi_npmm()}: \itemize{
#'  \item \strong{Known covariance} \itemize{
#'        \item \emph{diagonal covariance} We need the following additional
#'        arguments: \describe{
#'            \item{\code{cov_data}: a non-negative diagonal matrix, representing
#'                  the covariance of the data}{}
#'            \item{\code{prior_precision_scalar_eta}: a non-negative scalar,
#'                  representing the precision prior for the DP mean parameters}{}
#'            \item{\code{post_precision_scalar_eta}: initial value for the
#'                  posterior update of precision for the DP mean parameters}{}
#'            }
#'          \item \emph{full covariance} We need the following additional
#'        arguments: \describe{
#'            \item{\code{cov_data}: a positive definite matrix, representing
#'                  the covariance of the data}{}
#'            \item{\code{prior_cov_eta}: a positive definite matrix,
#'                  representing the covariance prior for the DP mean parameters}{}
#'            \item{\code{post_cov_eta}: initial value for the
#'                  posterior update of covariance for the DP mean parameters}{}
#'            }
#'          }
#'  \item \strong{Unknown covariance (Global)} \itemize{
#'    \item \emph{diagonal covariance} We need the following additional
#'    arguments: \describe{
#'      \item{\code{prior_shape_scalar_cov}: a non-negative scalar, representing
#'            the shape parameter of Gamma prior for the precision}{}
#'      \item{\code{prior_rate_scalar_cov}: a non-negative scalar, representing
#'            the rate parameter of Gamma prior for the precision}{}
#'      \item{\code{post_shape_scalar_cov}: initial value for posterior update of
#'            precision shape parameter}{}
#'      \item{\code{post_rate_scalar_cov}: initial value for posterior update of
#'            precision rate parameter}{}
#'      \item{\code{prior_precision_scalar_eta}: a non-negative scalar,
#'            representing the precision prior for the DP mean parameters}{}
#'      \item{\code{post_precision_scalar_eta}: initial value for the
#'            posterior update of precision for the DP mean parameters}{}
#'    }
#'    \item \emph{Inverse-Wishart} We need the following additional
#'    arguments: \describe{
#'      \item{\code{prior_df_cov}: a scalar as the degree of freedom parameter
#'            of the Inverse-Wishart prior, Default value D+2}{}
#'      \item{\code{prior_scale_cov}: positive-definite matrix as the scale
#'            parameter of the Inverse-Wishart prior}{}
#'      \item{\code{post_df_cov}: initial value for the posterior update of
#'            degree of freedom}{}
#'      \item{\code{post_scale_cov}: initial value for the posterior update of
#'            scale matrix}{}
#'      \item{\code{prior_cov_eta}: a positive definite matrix,
#'            representing the covariance prior for the DP mean parameters}{}
#'      \item{\code{post_cov_eta}: initial value for the
#'            posterior update of covariance for the DP mean parameters}{}
#'    }
#'    \item \emph{Cholesky-decomposition} We need the following additional
#'    arguments: \describe{
#'    \item{\code{prior_shape_diag_decomp}: a non-negative scalar as the shape
#'            parameter of Gamma prior for diagonal elements of the
#'            Cholesly-decomposed matrix}{}
#'    \item{\code{prior_rate_diag_decomp}: a non-negative scalar as the rate
#'            parameter of Gamma prior for diagonal elements of the
#'            Cholesly-decomposed matrix}{}
#'    \item{\code{prior_mean_offdiag_decomp}: a scalar as the mean
#'            parameter of Normal prior for off-diagonal elements of the
#'            Cholesly-decomposed matrix}{}
#'    \item{\code{prior_var_offdiag_decomp}: a non-negative scalar as the variance
#'            parameter of Normal prior for off-diagonal elements of the
#'            Cholesly-decomposed matrix}{}
#'    \item{\code{post_shape_diag_decomp}: initial value for posterior update
#'            of the shape parameter for diagonal elements}{}
#'    \item{\code{post_rate_diag_decomp}: initial value for posterior update
#'            of the rate parameter for diagonal elements}{}
#'    \item{\code{post_mean_offdiag_decomp}: initial value for posterior update
#'            of the mean parameter for off-diagonal elements}{}
#'    \item{\code{post_var_offdiag_decomp}: initial value for posterior update
#'            of the variance parameter for off-diagonal elements}{}
#'    \item{\code{prior_cov_eta}: a positive definite matrix,
#'            representing the covariance prior for the DP mean parameters}{}
#'    \item{\code{post_cov_eta}: initial value for the
#'            posterior update of covariance for the DP mean parameters}{}}
#'  }
#'  \item \strong{Unknown covariance (cluster-specific)} \itemize{
#'    \item \emph{Inverse Wishart} We need the following additional
#'    arguments: \describe{
#'      \item{\code{prior_df_cs_cov}: a vector representing degree of freedom
#'            parameters for each cluster-specific Inverse-Wishart prior}{}
#'      \item{\code{prior_scale_cs_cov}: an array of positive-definite matrices
#'            representing scale matrix parameters for each cluster-specific
#'            Inverse-Wishart prior}{}
#'      \item{\code{post_df_cs_cov}: initial value for posterior update of the
#'            degree of freedom parameters}{}
#'      \item{\code{post_scale_cs_cov}: initial value for posterior update of
#'            the scale matrix parameters}{}
#'      \item{\code{scaling_cov_eta}: a non-negative scaling factor for
#'            covariance matrix of the DP mean parameters}{}
#'    }
#'    \item \emph{Element-wise Gamma and Laplace prior} We need the following
#'    additional arguments: \describe{
#'      \item{\code{prior_shape_d_cs_cov}: a non-negative vector as shape
#'            parameters for cluster-specific Gamma priors of the diagonal
#'            elements}{}
#'      \item{\code{prior_rate_d_cs_cov}: a non-negative matrix as rate
#'            parameter for cluster-specific Gamma prior of the diagonal
#'            elements}{}
#'      \item{\code{prior_var_offd_cs_cov}: a non-negative vector as variance
#'            parameter for cluster-specific Laplace priors of the off-diagonal
#'            elements}{}
#'      \item{\code{post_shape_d_cs_cov}: initial value for posterior update of
#'            the diagonal shape parameters}{}
#'      \item{\code{post_rate_d_cs_cov}: initial value for posterior update of
#'            the diagonal rate parameters}{}
#'      \item{\code{post_var_offd_cs_cov}: initial value for posterior update of
#'            the off-diagonal variance parameters}{}
#'      \item{\code{scaling_cov_eta}: a non-negative scaling factor for
#'            covariance matrix of the DP mean parameters}{}
#'    }
#'    \item \emph{Element-wise Gamma and Normal prior} We need the following
#'    additional arguments: \describe{
#'      \item{\code{prior_shape_d_cs_cov}: a non-negative vector as shape
#'            parameters for cluster-specific Gamma priors of the diagonal
#'            elements}{}
#'      \item{\code{prior_rate_d_cs_cov}: a non-negative matrix as rate
#'            parameter for cluster-specific Gamma prior of the diagonal
#'            elements}{}
#'      \item{\code{prior_var_offd_cs_cov}: a non-negative scalar as variance
#'            parameter for cluster-specific Normal priors of the off-diagonal
#'            elements}{}
#'      \item{\code{post_shape_d_cs_cov}: initial value for posterior update of
#'            the diagonal shape parameters}{}
#'      \item{\code{post_rate_d_cs_cov}: initial value for posterior update of
#'            the diagonal rate parameters}{}
#'      \item{\code{post_mean_offd_cs_cov}: initial value for posterior update of
#'            the off-diagonal mean parameters}{}
#'      \item{\code{scaling_cov_eta}: a non-negative scaling factor for
#'            covariance matrix of the DP mean parameters}{}
#'    }
#'  }
#' }
#'
#' @param X input data as a matrix
#' @param variational_params number of clusters in the variational distribution
#' @param prior_shape_alpha shape parameter of Gamma prior for the DP
#' concentration parameter alpha. Default is 0.001
#' @param prior_rate_alpha rate parameter of Gamma prior for the DP
#' concentration parameter alpha. Default is 0.001
#' @param post_shape_alpha initial value for posterior update of shape parameter
#' for alpha. Default is 0.001
#' @param post_rate_alpha initial value for posterior update of ratee parameter
#' for alpha. Default is 0.001
#' @param prior_mean_eta mean vector of MVN prior for the DP mean parameters.
#' Default is zero vector
#' @param post_mean_eta initial value of posterior update for the DP mean
#' parameter
#' @param log_prob_matrix logarithm of cluster allocation probability matrix
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
#' @param ... additional parameters, further details given below
#'
#' @returns `[vimixr()]` returns a `list` with the following elements:
#'   - `alpha`: posterior DP concentration parameter
#'   - `Cluster number`: number of clusters from posterior probability allocation matrix
#'   - `Cluster Proportion`: cluster proportions from posterior probability allocation matrix
#'   - `log Probability matrix`: log of posterior probability allocation matrix
#'   - `ELBO`: Optimisation of the ELBO function
#'   - `Iterations`: Number of iterations required for convergence
#'   - `PCA_viz`: A PCA `[ggplot2]` plot to visualize the clustering of data based on cluster labels
#'   - `ELBO_viz`: A line `[ggplot2]` plot to visualize the ELBO optimisation
#'
#'
#' @importFrom Rfast rowsums colsums spdinv eachrow eachcol.apply Diag.fill
#' Diag.matrix
#' @importFrom ggplot2 ggplot aes geom_point geom_line labs
#' theme_minimal
#' @importFrom rlang .data
#' @importFrom stats prcomp
#'
#' @export
#'
#' @examples
#'
#' X <- rbind(matrix(rnorm(100, m=0, sd=0.5), ncol=2),
#'            matrix(rnorm(100, m=3, sd=0.5), ncol=2))
#'
#' #for fixed-diagonal
#' res <- cvi_npmm(X, variational_params = 20, prior_shape_alpha = 0.001,
#'          prior_rate_alpha = 0.001, post_shape_alpha = 0.001,
#'          post_rate_alpha = 0.001, prior_mean_eta = matrix(0, 1, ncol(X)),
#'          post_mean_eta = matrix(0.001, 20, ncol(X)),
#'          log_prob_matrix = t(apply(matrix(-3, nrow(X), 20), 1,
#'                              function(x){x/sum(x)})), maxit = 1000,
#'          fixed_variance = TRUE, covariance_type = "diagonal",
#'          prior_precision_scalar_eta = 0.001,
#'          post_precision_scalar_eta = matrix(0.001, 20, 1),
#'          cov_data = diag(ncol(X)))
#'  summary(res)
#'  plot(res)
#'

cvi_npmm <- function(X, variational_params,
                     prior_shape_alpha, prior_rate_alpha,
                     post_shape_alpha, post_rate_alpha,
                     prior_mean_eta, post_mean_eta,
                     log_prob_matrix,
                     maxit = 100,
                     covariance_type="full", fixed_variance=FALSE,
                     cluster_specific_covariance=TRUE,
                     variance_prior_type=c("IW", "decomposed", "sparse",
                                           "off-diagonal normal"),
                     ...
                     ){
  N <- nrow(X)
  D <- ncol(X)
  T0 <- variational_params
  log_prob_matrix <- log_prob_matrix[, order(Rfast::colsums(exp(log_prob_matrix)),
                                           decreasing = TRUE)]
  varargs <- list(...)
  params <- list()
  inverts <- list()

  params$N <- N #number of samples
  params$D <- D # number of variables
  params$T0 <- T0

  params$prior_mean_eta <- prior_mean_eta
  params$prior_shape_alpha <- prior_shape_alpha
  params$prior_rate_alpha <- prior_rate_alpha
  params$post_shape_alpha <- post_shape_alpha
  params$post_rate_alpha <- post_rate_alpha
  params$post_mean_eta <- post_mean_eta

  params$log_prob_matrix <- log_prob_matrix
  params$P <- t(apply(exp(log_prob_matrix), 1, function(x){x/sum(x)}))
  RP <- Rfast::colsums(params$P)


  #updating the parameter list based on the conditions
  if(covariance_type == "diagonal") {

    if(fixed_variance) {
      params$prior_precision_scalar_eta <- varargs$prior_precision_scalar_eta
      params$post_precision_scalar_eta <- varargs$post_precision_scalar_eta
      params$cov_data <- varargs$cov_data
      params_check(params, fixed_variance, covariance_type,
                   cluster_specific_covariance,
                   variance_prior_type)

      C00 <- diag(D)/varargs$prior_precision_scalar_eta #covariance of DP mean parameters
      inverts[["inv_C0"]] <- Rfast::spdinv(varargs$cov_data)
      inverts[["inv_C00"]] <- Rfast::spdinv(C00)

    } else {
      params$prior_shape_scalar_cov <- varargs$prior_shape_scalar_cov
      params$prior_rate_scalar_cov <- varargs$prior_rate_scalar_cov
      params$post_shape_scalar_cov <- varargs$post_shape_scalar_cov
      params$post_rate_scalar_cov <- varargs$post_rate_scalar_cov
      params$post_precision_scalar_eta <- varargs$post_precision_scalar_eta
      params$prior_precision_scalar_eta <- varargs$prior_precision_scalar_eta
      params_check(params, fixed_variance, covariance_type,
                   cluster_specific_covariance,
                   variance_prior_type)

      C00 <- diag(D)/varargs$prior_precision_scalar_eta #covariance of DP mean parameters
      inverts[["inv_C00"]] <- Rfast::spdinv(C00)

    }

  } else if(covariance_type == "full") {

    if(fixed_variance) {
      params$post_cov_eta <- varargs$post_cov_eta
      params$cov_data <- varargs$cov_data
      params$prior_cov_eta <- varargs$prior_cov_eta
      params_check(params, fixed_variance, covariance_type,
                   cluster_specific_covariance,
                   variance_prior_type)

      inverts[["inv_C0"]] <- Rfast::spdinv(varargs$cov_data)
      inverts[["inv_C00"]] <- Rfast::spdinv(varargs$prior_cov_eta)


    } else {
      if(!cluster_specific_covariance) {
        if(variance_prior_type == "IW"){
          params$prior_df_cov <- varargs$prior_df_cov
          params$prior_scale_cov <- varargs$prior_scale_cov
          params$post_df_cov <- varargs$post_df_cov
          params$post_scale_cov <- varargs$post_scale_cov
          params$post_cov_eta <- varargs$post_cov_eta
          params$prior_cov_eta <- varargs$prior_cov_eta
          params_check(params, fixed_variance, covariance_type,
                       cluster_specific_covariance,
                       variance_prior_type)

          inverts[["inv_V0"]] <- Rfast::spdinv(varargs$prior_scale_cov)
          inverts[["inv_C00"]] <- Rfast::spdinv(varargs$prior_cov_eta)


        } else if (variance_prior_type == "decomposed"){
          params$prior_shape_diag_decomp <- varargs$prior_shape_diag_decomp
          params$prior_rate_diag_decomp <- varargs$prior_rate_diag_decomp
          params$prior_mean_offdiag_decomp <- varargs$prior_mean_offdiag_decomp
          params$prior_var_offdiag_decomp <- varargs$prior_var_offdiag_decomp
          params$post_shape_diag_decomp <- varargs$post_shape_diag_decomp
          params$post_rate_diag_decomp <- varargs$post_rate_diag_decomp
          params$post_mean_offdiag_decomp <- varargs$post_mean_offdiag_decomp
          params$post_var_offdiag_decomp <- varargs$post_var_offdiag_decomp
          params$post_cov_eta <- varargs$post_cov_eta
          params$prior_cov_eta <- varargs$prior_cov_eta
          params_check(params, fixed_variance, covariance_type,
                       cluster_specific_covariance,
                       variance_prior_type)

          inverts[["inv_C00"]] <- Rfast::spdinv(varargs$prior_cov_eta)


        } else {
          stop("'variance_prior_type' can only be either 'IW' or 'decomposed'
               when 'cluster_specific_covariance' is FALSE")
        }

      }else{
        if(variance_prior_type == "IW"){
          params$prior_df_cs_cov <- varargs$prior_df_cs_cov
          params$prior_scale_cs_cov <- varargs$prior_scale_cs_cov
          params$post_df_cs_cov <- varargs$post_df_cs_cov
          params$post_scale_cs_cov <- varargs$post_scale_cs_cov
          params$scaling_cov_eta <- varargs$scaling_cov_eta
          params_check(params, fixed_variance, covariance_type,
                       cluster_specific_covariance,
                       variance_prior_type)


        } else if (variance_prior_type == "sparse"){
          params$prior_shape_d_cs_cov <- varargs$prior_shape_d_cs_cov
          params$prior_rate_d_cs_cov <- varargs$prior_rate_d_cs_cov
          params$prior_var_offd_cs_cov <- varargs$prior_var_offd_cs_cov
          params$post_shape_d_cs_cov <- varargs$post_shape_d_cs_cov
          params$post_rate_d_cs_cov <- varargs$post_rate_d_cs_cov
          params$post_var_offd_cs_cov <- varargs$post_var_offd_cs_cov
          params$scaling_cov_eta <- varargs$scaling_cov_eta
          params_check(params, fixed_variance, covariance_type,
                       cluster_specific_covariance,
                       variance_prior_type)


        } else if (variance_prior_type == "off-diagonal normal"){
          params$prior_shape_d_cs_cov <- varargs$prior_shape_d_cs_cov
          params$prior_rate_d_cs_cov <- varargs$prior_rate_d_cs_cov
          params$prior_var_offd_cs_cov <- varargs$prior_var_offd_cs_cov
          params$post_shape_d_cs_cov <- varargs$post_shape_d_cs_cov
          params$post_rate_d_cs_cov <- varargs$post_rate_d_cs_cov
          params$post_mean_offd_cs_cov <- varargs$post_mean_offd_cs_cov
          params$scaling_cov_eta <- varargs$scaling_cov_eta
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
  #plots
  pca <- prcomp(X)
  pca_df <- data.frame("PC1" = pca$x[,1],
                       "PC2" = pca$x[,2],
                       "Cluster" = as.factor(clustering))
  ggplot_pca <- ggplot2::ggplot(pca_df, ggplot2::aes(x = .data$PC1, y = .data$PC2,
                                                       color = .data$Cluster,
                                                     shape = .data$Cluster)) +
    ggplot2::geom_point(size = 3, alpha = 0.8) +
    ggplot2::labs(title = "PCA projection of clustered data", x = "PC 1", y = "PC 2") +
    ggplot2::theme_minimal()

  Elbo <- unlist(lapply(elbo_values[-1], sum))
  Elbo_df <- data.frame(x = 1:length(Elbo),
                        y = Elbo)
  ggplot_ELBO <- ggplot2::ggplot(Elbo_df, ggplot2::aes(x = .data$x,
                                                       y = .data$y)) +
    ggplot2::geom_line() +
    ggplot2::labs(title = "ELBO Optimisation", x = "Iterations", y = "ELBO") +
    ggplot2::theme_minimal()

  post_distribution = list()
  L1 <- params$post_mean_eta

  if(covariance_type == "diagonal") {

    if(fixed_variance) {
      L2 <- params$post_precision_scalar_eta

      post_distribution[["Mean"]] = L1/c(L2)

    } else {
      G1 <- params$post_shape_scalar_cov
      G2 <- params$post_rate_scalar_cov
      L2 <- params$post_precision_scalar_eta <- varargs$post_precision_scalar_eta

      post_distribution[["Mean"]] = L1/c(L2)
      post_distribution[["Precision"]] = (G1/G2)*diag(D)
      #G2/G1 instead of G1/G2 because prior on the precision diagonal scalar
    }

  } else if(covariance_type == "full") {

    if(fixed_variance) {
      L2 <- params$post_cov_eta
      L21 <- matrix(0, nrow = T0, ncol = D)
      for (i in 1:T0){
        L21[i,] = mat_mult(L1[i,, drop = FALSE], L2[,,i])
      }

      post_distribution[["Mean"]] = L21

    } else {
      if(!cluster_specific_covariance) {
        if(variance_prior_type == "IW"){
          nu <- params$post_df_cov
          V <- params$post_scale_cov
          L2 <- params$post_cov_eta

          L21 <- matrix(0, nrow = T0, ncol = D)
          for (i in 1:T0){
            L21[i,] = mat_mult(L1[i,, drop = FALSE], L2[,,i])
          }

          post_distribution[["Mean"]] = L21
          post_distribution[["Precision"]] = nu*V

        } else if (variance_prior_type == "decomposed"){
          a1 <- params$post_shape_diag_decomp
          b1 <- params$post_rate_diag_decomp
          mu1 <- params$post_mean_offdiag_decomp
          c1 <- params$post_var_offdiag_decomp
          L2 <- params$post_cov_eta

          mean_lower <- matrix(0, nrow = D, ncol = D) #mean matrix of the decomposed
          mean_lower[lower.tri(mean_lower, diag = FALSE)] <- mu1
          sigma_lower <- matrix(0, nrow = D, ncol = D) #var matrix of the decomposed
          sigma_lower[lower.tri(sigma_lower, diag = FALSE)] <- c1
          mean_L <- mean_lower + diag(sqrt(1/b1)*sqrt(pi)/beta(a1,0.5))
          diag(sigma_lower) <- (1/b1)*(a1 - (sqrt(pi)/beta(a1,0.5))^2)
          #expected inverse of C0; covariance matrix of data
          inv_C0 <- mat_mult(mean_L, t(mean_L)) +
            diag(Rfast::rowsums(sigma_lower))

          L21 <- matrix(0, nrow = T0, ncol = D)
          for (i in 1:T0){
            L21[i,] = mat_mult(L1[i,, drop = FALSE], L2[,,i])
          }

          post_distribution[["Mean"]] = L21
          post_distribution[["Precision"]] = inv_C0

        } else {
          stop("'variance_prior_type' can only be either 'IW' or 'decomposed'
               when 'cluster_specific_covariance' is FALSE")
        }

      }else{
        if(variance_prior_type == "IW"){
          nu1 <- params$post_df_cs_cov
          V1 <- params$post_scale_cs_cov

          V1_inv <- array(apply(V1, 3, function(x){spdinv(x)}), dim = dim(V1))
          #expectation of inverse of data covariance matrix
          inv_C0 <- sweep_3D(V1_inv, nu1, c(D, D, T0))

          post_distribution[["Mean"]] = L1
          post_distribution[["Precision_cs"]] = inv_C0

        } else if (variance_prior_type == "sparse"){
          a1 <- params$post_shape_d_cs_cov
          B1 <- params$post_rate_d_cs_cov

          #expectation of inverse of C0, data covariance matrix
          inv_C0 <- array(0, c(D, D, T0))
          for (i in 1:T0){
            inv_C0[,,i] <- temp <- Rfast::Diag.matrix(D, a1[1,i]/B1[i,])
          }

          post_distribution[["Mean"]] = L1
          post_distribution[["Precision_cs"]] = inv_C0
          post_distribution[["Precision_cs_a1"]] = a1
          post_distribution[["Precision_cs_B1"]] = B1

        } else if (variance_prior_type == "off-diagonal normal"){
          a1 <- params$post_shape_d_cs_cov
          B1 <- params$post_rate_d_cs_cov
          C1 <- params$post_mean_offd_cs_cov

          #expectation of inverse of data covariance matrix
          inv_C0 <- array(0, c(D, D, T0))
          for (i in 1:T0){
            inv_C0[,,i] <- temp <- Rfast::Diag.fill(C1[,,i], a1[1,i]/B1[i,])}

          post_distribution[["Mean"]] = L1
          post_distribution[["Precision_cs"]] = inv_C0

        } else {
          stop("'variance_prior_type' can only be either 'IW' or 'decomposed'
               when 'cluster_specific_covariance' is TRUE")
        }
      }

    }
  } else {
    stop("covariance_type can only be either 'diagonal' or 'full'.")
  }


  posterior <- c(list("alpha" = alpha0, "Cluster number" = clustnum,
                    "Cluster Proportion" = clust,
                    "log Probability matrix" = Plog), post_distribution)
  logBayes <- as.list(elbo_values[length(elbo_values)])$e_data - 
    as.list(elbo_values[1])$e_data
  optimisation <- list("ELBO" = elbo_values,
                       "Iterations" = (length(elbo_values)-1),
                       "logBF" = logBayes)

  output <-  list("posterior" = posterior, "optimisation" = optimisation,
                  "PCA_viz" = ggplot_pca,
                  "ELBO_viz" = ggplot_ELBO)
  class(output) <- "CVIoutput"

  return(output)
}

