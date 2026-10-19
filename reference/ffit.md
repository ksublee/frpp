# Fit a Flexible Residual Point Process

Estimate a composition model using an FRPP likelihood. Choose the
compatibility bridge or independent numeric parameter and optimizer
execution.

## Usage

``` r
ffit(object, inter_arrival, type = NULL, mark = NULL,
     N = NULL, Nc = NULL, N0 = NULL, Nc0 = NULL,
     lambda_component0 = "auto", estimate_residual = FALSE,
     optimizer = c("maxLik", "stats"), method = "BFGS",
     backend = c("bridge", "state"), reduced = TRUE,
     fixed = NULL, constraints = NULL, constraint = NULL,
     valid = NULL, ..., fit_backend = c("bridge", "independent", "auto"))
fitted_fspec(object)
finitial_state(object)
fdiagnostics(object)
```

## Arguments

- object:

  An `fspec` for fitting, or an `ffit` result for the helpers.

- inter_arrival, type, mark:

  Observed intervals, event types and optional marks, including an
  initial zero row.

- N, Nc:

  Optional cumulative count and mark-total paths.

- N0, Nc0:

  Initial counts and mark totals.

- lambda_component0:

  Initial kernel components. `"auto"` selects zero; specify a scalar or
  conformable matrix when known. Hawkes mean initialization is not used.

- estimate_residual:

  Keep the residual distribution fixed by default. Set `TRUE` to include
  its named parameters in the optimization.

- optimizer:

  Optional maxLik provider, or R's stats provider. The stats provider
  does not require maxLik.

- method:

  Optimization method, default BFGS. The stats provider also supports
  Nelder-Mead, CG and L-BFGS-B.

- backend:

  FRPP likelihood backend. `state` retains the numeric-kernel
  restrictions described in
  [`logLik.fspec`](https://ksublee.github.io/frpp/reference/operations.md).

- fit_backend:

  Named-only fitting execution choice, independent of `backend`. The
  default `bridge` uses public
  [`emhawkes::hfit`](https://ksublee.github.io/emhawkes/reference/hfit.html).
  `independent` registers numeric kernel and raw residual parameters
  with the public
  [`emhawkes::parameter_plan`](https://ksublee.github.io/emhawkes/reference/parameter_plan.html)
  API and calls maxLik or stats directly. It requires that development
  API and excludes event callbacks and legacy bound residual laws.
  Unsupported requests are errors. `auto` selects independent fitting
  when supported and otherwise uses the bridge, before optimization; it
  never retries a failed fit using another route.

- reduced:

  Share entries having equal starting values within each numeric kernel
  matrix. Function parameters sharing a name are shared under either
  setting.

- fixed:

  Logical mask or parameter positions, passed to the optimizer. With
  stats, names are also supported. The order follows emhawkes: mu,
  alpha, beta, eta, impact, mark density, residual density/CDF/quantile;
  duplicated function parameter names occur once. Fixed residual
  parameters are excluded entirely when `estimate_residual = FALSE`.

- constraints, constraint:

  Alternative names for optimization constraints. Supply at most one.
  Linear inequalities are `list(ineqA, ineqB)` with \\A\theta+B\ge0\\;
  columns follow the full parameter order, including fixed optimization
  parameters. Stats uses constrOptim and requires a strictly interior
  start for movable inequalities; it does not support equalities.

- valid:

  Optional function of the named full parameter vector, returning one
  nonmissing logical value. FALSE rejects a proposal with minus
  infinity. Fixed raw residual parameters are also included. Validator
  and callback errors propagate. Supply distribution-specific domain
  checks when needed.

- ...:

  Named optimizer arguments, including `control`, `grad`, `hess`, and
  stats box bounds. For stats use `control$maxit`; maxLik uses
  `control$iterlim`. FRPP objective, analytic Hawkes gradient and Hawkes
  final-Hessian selection cannot be overridden through these arguments.

## Details

The bridge delegates fitting execution to public emhawkes functions.
Independent fitting owns observation preparation, optimization and
fitted-model assembly; it shares only the public parameter contract and
the selected likelihood's public kernel/residual APIs. Its stats adapter
retains the bridge's fixed/free, constraint reduction, finite-difference
and covariance rules. Neither route uses internal emhawkes functions or
Hawkes curvature stabilization. The result's `fitting` field records
requested and used fitting routes and an automatic fallback reason.
Earlier saved results may lack this field. Likelihood evaluation and
default BFGS fitting do not generate random values. Stochastic
optimizers can consume randomness. User callbacks must be deterministic.

Residual parameters are fixed by default to avoid silently estimating
two confounded scales. Joint estimation requires an identifiable model:
fixing or normalizing the relevant scale remains the user's modeling
choice. This API does not automatically diagnose all forms of
nonidentification or normalize arbitrary laws. Matrix sharing is
preserved on reconstruction; residual parameters shared among density,
CDF, quantile and generator functions are updated together.

The result supports coef, logLik, nobs, vcov, summary and fitted.
`fitted_fspec` and `fitted` return the estimated composition model;
`finitial_state` returns the initialization arguments for likelihood
replay or simulation. Save the result with saveRDS and restore it with
readRDS. Captured callback environments remain references; external
state is not frozen.

Wald covariance is withheld on optimizer failure, missing or nonpositive
curvature, near detected constraint boundaries, and equality
constraints. The provider result remains in `fit$fit`: maxLik retains
its class; the independent stats result has class `frpp_stats_fit`,
supporting coef, logLik, vcov, summary and print. Bridge stats results
retain their existing class. A positive covariance is reported as
`available_unvalidated`: numerical step sensitivity, nonsmooth
densities, weak identification and actual confidence-interval coverage
have not been validated. A supplied analytic Hessian can avoid
finite-difference error for a smooth model. No eigenvalue clipping or
absolute-variance correction is used. The tzexp transition may give a
nonsmooth optimum; consider profile likelihood or bootstrap separately
rather than interpreting a success code as valid inference.

## Value

`ffit` returns a list of class `ffit`, including the provider fit,
estimated fspec, coefficients, covariance, likelihood, initialization
and diagnostics. `fitted_fspec` returns an fspec, `finitial_state` a
list, and `fdiagnostics` a list containing convergence and
Wald-availability information. Likelihood degrees of freedom count only
active optimization parameters; nobs excludes the initial row.
Unavailable free-parameter covariance entries are NA.

## Examples

``` r
model <- fspec(methods::new("hspec", mu = 1, alpha = 0, beta = 1),
  residual_law(c(rate = 1), stats::dexp, stats::pexp, stats::qexp, stats::rexp))
fit <- ffit(model, c(0, .2, .3, .5), optimizer = "stats",
  fixed = c(FALSE, TRUE, TRUE), lambda_component0 = 0)
coef(fit)
#>    mu1 alpha1  beta1 
#>      3      0      1 
summary(fit)
#> FRPP maximum likelihood fit (stats, BFGS, bridge)
#> Log likelihood: 0.2958369
#> Convergence code: 0 (Normal convergence.)
#>        Estimate Std. Error z value Pr(>|z|)  
#> mu1      3.0000     1.7321  1.7321  0.08326 .
#> alpha1   0.0000     0.0000      NA       NA  
#> beta1    1.0000     0.0000      NA       NA  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> Fixed parameters: alpha1, beta1
#> Wald covariance: available_unvalidated
fdiagnostics(fit)
#> $optimizer
#> [1] "stats"
#> 
#> $method
#> [1] "BFGS"
#> 
#> $backend
#> [1] "bridge"
#> 
#> $optimizer_code
#> [1] 0
#> 
#> $optimizer_message
#> [1] "Normal convergence."
#> 
#> $optimizer_success
#> [1] TRUE
#> 
#> $estimate_residual
#> [1] FALSE
#> 
#> $wald
#> [1] "available_unvalidated"
#> 
#> $reason
#> NULL
#> 
#> $active
#> [1]  TRUE FALSE FALSE
#> 
#> $note
#> [1] "Optimizer curvature is not a validation of FRPP Wald coverage; no Hawkes Hessian stabilization is applied."
#> 
estimated <- fitted(fit)
do.call(logLik, c(list(object = estimated, inter_arrival = c(0, .2, .3, .5)),
  finitial_state(fit)))
#> loglikelihood 
#>     0.2958369 
```
