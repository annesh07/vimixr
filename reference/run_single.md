# CVI implementation for one set of initial parameters

CVI implementation for one set of initial parameters

## Usage

``` r
run_single(
  config,
  X,
  N,
  D,
  T0,
  prior_shape_alpha,
  prior_rate_alpha,
  post_shape_alpha,
  post_rate_alpha,
  prior_mean_eta,
  post_mean_eta,
  fixed_variance,
  covariance_type,
  cluster_specific_covariance,
  variance_prior_type,
  maxit,
  varargs
)
```

## Arguments

- config:

  List of inputs that are generated if not user-provided

- X:

  the data matrix

- N:

  samples of X

- D:

  dimensions of X

- T0:

  variational clusters

- prior_shape_alpha:

  shape parameter of Gamma prior for the DP concentration parameter
  alpha. Default is 0.001

- prior_rate_alpha:

  rate parameter of Gamma prior for the DP concentration parameter
  alpha. Default is 0.001

- post_shape_alpha:

  initial value for posterior update of shape parameter for alpha.
  Default is 0.001

- post_rate_alpha:

  initial value for posterior update of ratee parameter for alpha.
  Default is 0.001

- prior_mean_eta:

  mean vector of MVN prior for the DP mean parameters. Default is zero
  vector

- post_mean_eta:

  initial value of posterior update for the DP mean parameter

- fixed_variance:

  covariance matrix of the data is considered known (fixed) or unknown.

- covariance_type:

  covariance matrix is considered diagonal or full.

- cluster_specific_covariance:

  covariance matrix is specific to a cluster allocation or it is same
  over all cluster choices.

- variance_prior_type:

  For unknown and full covariance matrix, choice of matrix prior is
  either Inverse-Wishart ('IW') or Cholesky-decomposed ('decomposed').
  For unknown, full and cluster-specific covariance matrix, choice of
  matrix prior is either Inverse-Wishart ('IW'), element-wise Gamma and
  Laplace distributed ('sparse') or element-wise Gamma and Normal
  distributed ('off-diagonal normal')

- maxit:

  Maximum number of iterations for variational updates

- varargs:

  List of case specific parameters

## Value

a `list` with the following elements:

- `alpha`: posterior DP concentration parameter

- `Cluster number`: number of clusters from posterior probability
  allocation matrix

- `Cluster Proportion`: cluster proportions from posterior probability
  allocation matrix

- `log Probability matrix`: log of posterior probability allocation
  matrix

- `ELBO`: Optimisation of the ELBO function

- `Iterations`: Number of iterations required for convergence
