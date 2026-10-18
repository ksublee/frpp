# Bundle Residual Distribution Functions

Store a named parameter vector and residual distribution functions
without evaluating the functions or drawing random numbers.

## Usage

``` r
residual_law(param = numeric(), density = NULL, cdf = NULL,
             quantile = NULL, rng = NULL, name = "custom")
```

## Arguments

- param:

  Finite numeric vector with unique, nonempty parameter names.

- density, cdf, quantile, rng:

  Distribution functions, or `NULL`. Use R distribution signatures: the
  variate, probability or sample count is the first argument. Named
  parameters are matched to explicit formal arguments as in
  [`emhawkes::set_residual`](https://ksublee.github.io/emhawkes/reference/set_residual.html).
  Log-density and tail options are supported there.

- name:

  A descriptive label, not a verified distribution-family identifier.

## Details

At least a density or random generator is required. Simulation needs a
generator; likelihood evaluation needs a density, and multivariate
likelihoods also need a CDF. Residual inference needs a CDF and quantile
function. The constructor does not establish consistency of
user-supplied functions, normalization of residual scale, existence of
moments, or parameter identification.

## Value

An object of class `frpp_residual_law`.
