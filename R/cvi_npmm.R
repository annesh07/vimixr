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
#' @param log_prob_matrix logarithm of cluster allocation probability matrix. 
#' Default is NULL
#' @param maxit maximum number of iterations. Default is 100
#' @param n_inits Number of random initialisations if log_prob_matrix and other 
#' case-specific hyperparameters are NULL. Default is 5
#' @param Seed Seeds for random initialisation; either a vector of n_inits 
#' integers or NULL. Default is NULL.
#' @param parallel Logical input for parallelisation. Default is FALSE
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
#' @importFrom Rfast rowsums colsums spdinv eachrow eachcol.apply Diag.fill  Diag.matrix
#' @importFrom ggplot2 ggplot aes geom_point geom_line labs theme_minimal
#' @importFrom rlang .data
#' @importFrom irlba prcomp_irlba
#' @importFrom parallel detectCores makeCluster stopCluster clusterExport
#' @importFrom utils tail
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
#'                              function(x){x/sum(x)})), maxit = 100,
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
                     log_prob_matrix = NULL,
                     maxit = 100,
                     n_inits = 5,
                     Seed = NULL,
                     parallel = FALSE,
                     covariance_type="full", fixed_variance=FALSE,
                     cluster_specific_covariance=TRUE,
                     variance_prior_type=c("IW", "decomposed", "sparse",
                                           "off-diagonal normal"),
                     ...
                     ){
  N <- nrow(X)
  D <- ncol(X)
  T0 <- variational_params
  varargs <- list(...)
  
  #if multi-initialisation is required
  need_random_logP <- is.null(log_prob_matrix)
  sparse_case <- (covariance_type == "full" && !fixed_variance && cluster_specific_covariance && variance_prior_type == "sparse")
  need_random_priors <- sparse_case && (is.null(varargs$prior_shape_d_cs_cov) || is.null(varargs$prior_rate_d_cs_cov))
  effective_n_inits <- if (need_random_logP || need_random_priors) n_inits else 1L
  
  #vector of seeds if random generalisation required
  if (need_random_logP || need_random_priors){
    if (length(Seed)==0) Seed <- sample.int(1e+7, n_inits, replace = FALSE)
  } else {
    Seed = "No random initialisation used"
  }
  
  # Unified list of inputs as configs
  empBayes_values <- rep(0, effective_n_inits)
  configs <- vector("list", effective_n_inits)
  for (i in 1:effective_n_inits){
    seed0 <- Seed[i]
    config_i <- list()
    if (need_random_logP) {
      logP <- generate_log_prob(N, T0, seed0)
    } else {
      logP <- log_prob_matrix
    }
    config_i$log_prob_matrix <- logP
    if (need_random_priors) {
      cs_priors <- eBa0(logP, X)
      empBayes_values[i] <- cs_priors
      config_i$prior_shape_d_cs_cov <- matrix(cs_priors, 1, T0)
      config_i$prior_rate_d_cs_cov <- matrix(cs_priors, T0, D)
    } else {
      config_i$prior_shape_d_cs_cov <- varargs$prior_shape_d_cs_cov
      config_i$prior_rate_d_cs_cov <- varargs$prior_rate_d_cs_cov
    }
    
    configs[[i]] <- config_i
  }
  
  #wrapping the cvi function
  cvi_wrapper <- function(config) {
    run_single(config, 
               X = X, N = N, D = D, T0 = T0,
               prior_shape_alpha = prior_shape_alpha, 
               prior_rate_alpha = prior_rate_alpha,
               post_shape_alpha = post_shape_alpha, 
               post_rate_alpha = post_rate_alpha,
               prior_mean_eta = prior_mean_eta, 
               post_mean_eta = post_mean_eta,
               fixed_variance = fixed_variance, 
               covariance_type = covariance_type,
               cluster_specific_covariance = cluster_specific_covariance,
               variance_prior_type = variance_prior_type, 
               maxit = maxit, 
               varargs = varargs)
  }
  # Parallel setup if requested and multi-init
  if (effective_n_inits > 1 && parallel &&
      requireNamespace("pbapply", quietly = TRUE)) {
    
    n_cores <- parallel::detectCores()
    if (.Platform$OS.type == "unix") {
      results <- pbapply::pblapply(configs, cvi_wrapper, cl = n_cores)
      
    } else {
      n_cores <- max(1, n_cores - 1)
      cl <- parallel::makeCluster(n_cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterExport(cl, varlist = ls(environment()), envir = environment())
      results <- pbapply::pblapply(configs, cvi_wrapper, cl = cl)
    }
    
  } else {
    results <- lapply(configs, cvi_wrapper)
  }
  

  # Select best based on VLL if multi-init; else just the single result
  if (effective_n_inits > 1){
    final_elbos <- sapply(results, 
                          function(out){as.numeric(utils::tail(out$optimisation$ELBO, 1)[[1]]["e_data"])})
    best_idx <- which.max(final_elbos)
    best_result <- results[[best_idx]]
    best_result$index <- best_idx
  } else {
    best_result <- results[[1]]
  }
  # Extract data for plots
  posterior <- best_result$posterior
  optimisation <- best_result$optimisation
  clustering <- apply(posterior[["log Probability matrix"]], 1, which.max)
  
  # PCA plot
  pca <- irlba::prcomp_irlba(X,2)
  #variation explained
  var_explained <- pca$sdev^2 / pca$totalvar
  pc1_pct <- round(var_explained[1] * 100, 2)
  pc2_pct <- round(var_explained[2] * 100, 2)
  #the plot
  pca_df <- data.frame("PC1" = pca$x[,1], "PC2" = pca$x[,2], "Cluster" = as.factor(clustering))
  ggplot_pca <- ggplot2::ggplot(pca_df, ggplot2::aes(x = .data$PC1, y = .data$PC2, color = .data$Cluster, shape = .data$Cluster)) +
    ggplot2::geom_point(size = 3, alpha = 0.8) +
    ggplot2::labs(title = "PCA projection of SparseDPMM clusters", 
                  x = paste0("PC 1 (", pc1_pct, "%)"), 
                  y = paste0("PC 2 (", pc2_pct, "%)")) +
    ggplot2::theme_minimal()
  
  # ELBO plot
  Elbo <- unlist(lapply(optimisation$ELBO[-1], sum))
  Elbo_df <- data.frame(x = 1:length(Elbo), y = Elbo)
  ggplot_ELBO <- ggplot2::ggplot(Elbo_df, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_line() +
    ggplot2::labs(title = "ELBO Optimisation", x = "Iterations", y = "ELBO") +
    ggplot2::theme_minimal()
  
  # Attach to best_result
  best_result$PCA_viz <- ggplot_pca
  best_result$ELBO_viz <- ggplot_ELBO
  
  #Attach seeds used for reproducibility
  best_result$Seed_used <- Seed
  
  #Attach the empirical Bayes based output for a0 hyper-parameter of Sparse DPMM
  if (need_random_priors) {
    best_result$Empirical_Bayes_estimates <- empBayes_values
  }
  
  class(best_result) <- "CVIoutput" 
  return(best_result)  
}

