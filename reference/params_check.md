# Function to check the list of type-specific arguments

Function to check the list of type-specific arguments

## Usage

``` r
params_check(
  params,
  fixed_variance = FALSE,
  covariance_type = "diagonal",
  cluster_specific_covariance = TRUE,
  variance_prior_type = c("IW", "decomposed", "sparse", "off-diagonal normal")
)
```

## Arguments

- params:

  the list of required parameters

- fixed_variance:

  whether covariance is assumed fixed or not; can be TRUE or FALSE

- covariance_type:

  structure of covariance matrix; can be "diagonal" or "full"

- cluster_specific_covariance:

  whether covariance matrix is cluster specific or not; can be TRUE or
  FALSE

- variance_prior_type:

  prior distribution for the covariance matrix; can be "IW" or
  "decomposed" when cluster_specific_covariance = FALSE, or can be "IW",
  "sparse" or "off-diagonal normal" otherwise

## Value

stops the code if the required list of arguments are not present
