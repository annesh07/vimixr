# Root for a0 hyper-parameter for Sparse DPMM

Root for a0 hyper-parameter for Sparse DPMM

## Usage

``` r
eBa0(
  logP,
  X,
  a_min = min(1e-08, 1/ncol(X)),
  a_max = max(1e+06, ncol(X)),
  grid_points = min(ncol(X), 10000)
)
```

## Arguments

- logP:

  log of probability allocation matrix

- X:

  observed data

- a_min:

  minimum value of a0 for grid search

- a_max:

  maximum value of a0 for grid search

- grid_points:

  number of points for grid search

## Value

No return value, called for side effects.
