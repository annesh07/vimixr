# Update of the variational parameters

Update of the variational parameters

## Usage

``` r
CVI_update_function(
  fixed_variance = FALSE,
  covariance_type = "diagonal",
  cluster_specific_covariance = TRUE,
  variance_prior_type = c("IW", "decomposed", "sparse", "off-diagonal normal"),
  X,
  inverts,
  params
)
```

## Arguments

- fixed_variance:

  whether the covariance is fixed or estimated. Default is `FALSE` which
  means it is estimated.

- covariance_type:

  The assumed type of the covariance matrix. Can be either `"diagonal"`
  if it is the identify multiplied by a scalar, or `"full"` for a fully
  unspecified covariance matrix.

- cluster_specific_covariance:

  whether the the covariance is shared across estimated clusters or is
  cluster specific. Default is `TRUE` which means it is cluster
  specific.

- variance_prior_type:

  character string specifying the type of prior distribution for the
  covariance when cluster_specific_covariance is `TRUE`. Can be either
  `"IW"` or `"decomposed"` if `cluster_specific_covariance` is `FALSE`,
  and can be either `"IW"`, `"sparse"` or `"off-diagonal normal"`
  otherwise.

- X:

  the data matrix

- inverts:

  a list of inverses

- params:

  a list of required arguments

## Value

Updated parameters
