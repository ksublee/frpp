# Get started

## Installation

Install the development package from GitHub:

``` r

install.packages("remotes")
remotes::install_github("ksublee/frpp", build_vignettes = FALSE)
```

The default backend requires emhawkes 1.0.0 or later. The optional
`state` backend uses the public shared-kernel API in GitHub development
emhawkes 1.0.1.9000. Install that development version when using state
calculations:

``` r

remotes::install_github("ksublee/emhawkes", build_vignettes = FALSE)
```

The default bridge works without that API. State requests report an
explicit error when the API is absent. See the [shared-kernel
contract](https://ksublee.github.io/emhawkes/articles/shared-kernel.html).

## Define a model

A model combines an emhawkes kernel with a residual law. The following
law has Gamma shape 2 and rate 2, hence mean 1:

``` r

kernel <- methods::new("hspec", mu = 0.3, alpha = 0.2, beta = 1)
law <- residual_law(c(shape = 2, rate = 2),
  density = stats::dgamma, cdf = stats::pgamma,
  quantile = stats::qgamma, rng = stats::rgamma, name = "Gamma(2,2)")
model <- fspec(kernel, law)
model
#> FRPP specification (S4 fspec)
#> Event types: 1
#> Residual law: Gamma(2,2)
#> Residual candidates refreshed after every global event.
```

The kernel supplies the baseline and exponential state dynamics. The
residual law supplies the density, distribution, quantile and
random-generation functions. Constructing these objects does not run the
residual callbacks.

For simulation, the current implementation requires equal decay rates
within each row of the kernel matrix. Numeric marked effects and kernel
mappings are supported; more general callbacks use the compatibility
backend.

## Simulate and inspect a path

``` r

path <- fsim(model, size = 1000, seed = 73, lambda_component0 = 0)
head(data.frame(time = cumsum(path$inter_arrival),
  type = path$type, mark = path$mark))
#>       time type mark
#> 1 0.000000    0    0
#> 2 2.212690    1    1
#> 3 4.776757    1    1
#> 4 5.266585    1    1
#> 5 5.848559    1    1
#> 6 8.591268    1    1
```

`size` includes the initial row, so this request generates 999 events.
`inter_arrival`, `type` and `mark` start with zero. Event types index
rows of the driver, while `type_col_map` selects kernel columns updated
after each event. `N` and `Nc` contain event counts and accumulated
marks.

`lambda_component0` contains initial intensity components **excluding
the baseline**. A scalar expands to every component; a matrix must match
the kernel shape. FRPP’s automatic rule selects zero, with a message
recommending an explicit state. A Hawkes stationary mean is not assumed
for a general residual law.

The standard interface returns the same path:

``` r

same_path <- stats::simulate(model, size = 1000, seed = 73,
  lambda_component0 = 0)
stopifnot(identical(path, same_path))
```

With an explicit `seed`, the caller’s RNG kind and state are restored,
including on errors. With `seed = NULL`, simulation advances the
caller’s RNG normally.

## Evaluate the likelihood

``` r

logLik(model, path$inter_arrival, path$type, mark = path$mark,
  lambda_component0 = 0)
#> loglikelihood 
#>     -1909.699
```

Include the initial row when supplying observations. For a univariate
model, `type` may be omitted. Use the same initial components and
initial counts used for the observation window.

With `infer = TRUE`, the result also contains driver states, component
histories and residual diagnostics. These driving quantities are not
automatically the conditional intensities of a general FRPP.

Continue with [Fitting and
inference](https://ksublee.github.io/frpp/articles/fitting-and-inference.md),
or see the [multivariate
examples](https://ksublee.github.io/frpp/articles/frpp-examples.md).

## Optional independent simulator

The default simulator remains the compatibility bridge. With current
GitHub development emhawkes, choose
`fsim(model, backend = "state", ...)` to use the independent numeric
FRPP event loop. It preserves global residual refresh and path/RNG
behavior, but is currently slower on larger paths. Event callbacks
remain in the bridge. `backend = "auto"` opts into capability-based
selection; see [Model and
computation](https://ksublee.github.io/frpp/articles/model-and-computation.md)
for the support contract.
