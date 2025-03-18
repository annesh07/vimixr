test_that("Number of clusters from latent allocation:",{

  expect_equal(cvi_npmm(X = rbind(matrix(rnorm(100, m=0, sd=0.5), ncol=2),
                                  matrix(rnorm(100, m=3, sd=0.5), ncol=2)),
                        variational_params = 20, prior_shape_alpha = 0.001,
                        prior_rate_alpha = 0.001, post_shape_alpha = 0.001,
                        post_rate_alpha = 0.001, prior_mean_eta = matrix(0, 1, 2),
                        post_mean_eta = matrix(0.001, 20, 2),
                        log_prob_matrix = t(apply(matrix(0.001, 100, 20), 1,
                                                  function(x){x/sum(x)})),
                        maxit = 1000,
                        covariance_type="diagonal",fixed_variance=TRUE,
                        prior_precision_scalar_eta = 0.001,
                        post_precision_scalar_eta = matrix(0.001, 20, 1),
                        cov_data = diag(2))$posterior$'Cluster number', 2)})
