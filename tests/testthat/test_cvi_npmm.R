test_that("Number of clusters from latent allocation:",{
  X <- rbind(matrix(rnorm(100, m=0, sd=0.5), ncol=2),
             matrix(rnorm(100, m=3, sd=0.5), ncol=2))
  prior_precision_scalar_eta <- 0.001
  post_precision_scalar_eta <- matrix(0.001, 20, 1)
  cov_data <- diag(ncol(X))
  expect_equal(cvi_npmm(X, variational_params = 20, prior_shape_alpha = 0.001,
                        prior_rate_alpha = 0.001, post_shape_alpha = 0.001,
                        post_rate_alpha = 0.001, prior_mean_eta = matrix(0, 1, ncol(X)),
                        post_mean_eta = matrix(0.001, 20, ncol(X)),
                        log_prob_matrix = t(apply(matrix(0.001, nrow(X), 20), 1,
                                                  function(x){x/sum(x)})), maxit = 1000,
                        covariance_type="diagonal",fixed_variance=TRUE,
                        prior_precision_scalar_eta,
                        post_precision_scalar_eta,
                        cov_data)$posterior$'Cluster number', 2)})
