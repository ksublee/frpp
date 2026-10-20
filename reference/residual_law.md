# Bundle Residual Distribution Functions

Store a named parameter vector and residual distribution functions
without evaluating the functions or drawing random numbers.

## Usage

``` r
residual_law(param = numeric(), dresidual = NULL, presidual = NULL,
             qresidual = NULL, rresidual = NULL, name = "custom",
             density = NULL, cdf = NULL, quantile = NULL, rng = NULL)
```

## Arguments

- param:

  Finite numeric vector with unique, nonempty parameter names.

- dresidual, presidual, qresidual, rresidual:

  Residual density, cumulative distribution, quantile and
  random-generation functions, or `NULL`. Use R distribution signatures:
  the variate, probability or sample count is the first argument. Named
  parameters are matched to explicit formal arguments as in
  [`emhawkes::set_residual`](https://ksublee.github.io/emhawkes/reference/set_residual.html).
  Log-density and tail options are supported there.

- density, cdf, quantile, rng:

  Compatibility aliases for `dresidual`, `presidual`, `qresidual` and
  `rresidual`, respectively. Supply at most one non-NULL function per
  pair; specifying both is an error, even when they are the same
  function. No deprecation warning is emitted.

- name:

  A descriptive label, not a verified distribution-family identifier.

## Details

Supply only the functions needed for your operation. Construction
requires at least `dresidual` or `rresidual`. Simulation needs only
`rresidual`; univariate likelihood and fitting need only `dresidual`.
Multivariate likelihood and fitting additionally need `presidual`.
`infer = TRUE` needs `dresidual`, `presidual` and `qresidual`, and does
not require `rresidual`. Missing functions are reported when an
operation requires them, before callbacks or random draws. The original
positional order of the first six arguments is unchanged. Stored law
objects keep their existing representation, and callback functions,
environments and defaults are preserved. The constructor does not
establish consistency of user-supplied functions, normalization of
residual scale, existence of moments, or parameter identification.

## Value

An object of class `frpp_residual_law`.

## Examples

``` r
# Simulation needs only a residual generator.
law <- residual_law(c(rate = 1), rresidual = stats::rexp)
model <- fspec(methods::new("hspec", mu = 1, alpha = 0, beta = 1), law)
fsim(model, size = 4, seed = 73, lambda_component0 = 0)
#> FRPP path: 3 events, 1 event types.

# Density alone is enough for a univariate likelihood.
law <- residual_law(c(rate = 1), dresidual = stats::dexp)
model <- fspec(methods::new("hspec", mu = 1, alpha = 0, beta = 1), law)
stats::logLik(model, c(0, .2, .3), lambda_component0 = 0)
#> loglikelihood 
#>          -0.5 
```
