# Simulate and Evaluate the Composition Prototype

Use a compatibility bridge, or assemble a numeric FRPP likelihood from
public emhawkes kernel replay output and the residual law.

## Usage

``` r
fsim(object, size = 100, lambda_component0 = "auto", N0 = NULL,
    Nc0 = NULL, seed = NULL, ...)
# S3 method for class 'fspec'
simulate(object, nsim = 1, seed = NULL, size = 100,
    lambda_component0 = "auto", N0 = NULL, Nc0 = NULL, ...)
# S3 method for class 'fspec'
logLik(object, inter_arrival, type = NULL, mark = NULL,
    N = NULL, Nc = NULL, N0 = NULL, Nc0 = NULL,
    lambda_component0 = "auto", infer = FALSE,
    backend = c("bridge", "state"), ...)
```

## Arguments

- object:

  An `fspec`. Old `frpp_spec` objects are also supported.

- nsim:

  One path is supported. Use `size` for its length.

- seed:

  Optional seed. When supplied, the caller's RNG state is restored.

- size:

  Number of output rows including the initial row.

- lambda_component0:

  Initial kernel components. Automatic initialization uses zero for the
  residual-distribution model. Supply an explicit state when known. A
  Hawkes stationary mean is not automatically applied to FRPP.

- N0, Nc0:

  Initial counts and mark totals, passed to emhawkes.

- inter_arrival, type, mark:

  Observations starting with an initial zero row.

- N, Nc:

  Optional observed count paths.

- infer:

  Return driver states, likelihood and residual diagnostics.

- backend:

  `bridge` delegates to the existing FRPP likelihood. `state`
  independently assembles the residual likelihood from the kernel's
  public replay output. This requires the emhawkes development kernel
  API; with earlier versions, use `bridge`.

- ...:

  Additional public emhawkes arguments for the bridge. The state
  prototype accepts none.

## Details

Simulation retains the existing R FRPP loop and requires equal decay
rates within each row. Every global event refreshes all candidate
residuals. The state backend currently supports numeric kernels, numeric
mark effects and kernel mappings, with equal decay rates within each row
and no impact or mark-density callbacks. Callback specifications remain
available through the bridge. Kernel states and their integrals are
drivers; with a non-exponential residual law they are not automatically
the FRPP conditional intensities.

For event type j in a global inter-event interval, its contribution is
\$\$\log v_j + \log f(A_j) + \sum\_{k \ne j}\log S(A_k),\$\$ where v is
the pre-event kernel driver and A is its interval integral. The public
replay call returns kernel-only interval integrals and pre-event sums.
The FRPP layer adds baseline contributions. Full component histories are
retained only with `infer = TRUE`; interval integrals and sums use
memory proportional to path length times the number of event types. This
backend does not evaluate a Hawkes likelihood. `fsim()` is the primary
simulation function; the
[`simulate()`](https://rdrr.io/r/stats/simulate.html) methods call it
and return the same path.

The existing `emhawkes.engine` option retains its original meaning. In
particular, the bridge does not support a forced C++ FRPP likelihood.

## Value

Simulation returns a `frpp_path` list containing its composition model
and the original path fields other than the transitional `hspec` field.
Likelihood evaluation returns the existing named numeric value, or an
inference list. The state backend additionally returns
`integrated_driver`. No fitted degrees of freedom or confidence
intervals are supplied.

## Examples

``` r
kernel <- methods::new("hspec", mu = .3, alpha = .2, beta = 1)
model <- fspec(kernel, residual_law(c(shape = 2, rate = 2),
    stats::dgamma, stats::pgamma, stats::qgamma, stats::rgamma))
path <- fsim(model, size = 12, seed = 73, lambda_component0 = 0)
stats::logLik(model, path$inter_arrival, path$type, lambda_component0 = 0)
#> loglikelihood 
#>     -22.68649 
if ("exp_kernel_replay" %in% getNamespaceExports("emhawkes")) {
  stats::logLik(model, path$inter_arrival, path$type, lambda_component0 = 0,
                backend = "state")
}
```
