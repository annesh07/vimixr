# sparse_cov_op

Calculate the sum, squared sum and log sum of off-diagonal vector
elements from the covariance array

## Usage

``` r
sparse_cov_op(X, P, inv_C0, L1)
```

## Arguments

- X:

  data matrix

- P:

  probability allocation matrix

- inv_C0:

  matrix corresponding diagonal elements of the cluster precision
  matrices

- L1:

  cluster mean matrix

## Value

likelihood term calculation in elbo
