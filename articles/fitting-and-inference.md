# Fitting and inference

## A small identifiable example

Unit-exponential residuals recover a Hawkes likelihood. With zero
excitation, only the baseline is identifiable: the decay rate has no
effect. Fix the excitation and decay parameters in this example.

``` r

model <- fspec(
  methods::new("hspec", mu = 1, alpha = 0, beta = 1),
  residual_law(c(rate = 1), stats::dexp, stats::pexp,
    stats::qexp, stats::rexp))
dt <- c(0, 0.2, 0.3, 0.5)
fit <- ffit(model, dt, optimizer = "stats", lambda_component0 = 0,
  fixed = c(FALSE, TRUE, TRUE), control = list(reltol = 1e-12))
coef(fit)
#>    mu1 alpha1  beta1 
#>      3      0      1
```

There are three events in one unit of time; the baseline MLE is 3. The
tighter optimizer tolerance makes this example easier to compare with
that analytic answer; it does not validate inference for more
complicated models.

## Optimizers and parameter sharing

The default `optimizer = "maxLik"` requires the optional maxLik package.
`optimizer = "stats"` uses R’s
[`optim()`](https://rdrr.io/r/stats/optim.html) or
[`constrOptim()`](https://rdrr.io/r/stats/constrOptim.html) through
emhawkes, without requiring maxLik. Controls differ: stats uses
`control$maxit`, while maxLik uses `control$iterlim`.

By default, equal starting values within a numeric matrix share a
parameter. Use `reduced = FALSE` for separate entries. Function
parameters with the same name remain shared. Inspect `coef(fit)` to see
the resulting names and order before specifying fixed parameters or
constraints.

`constraint` and `constraints` are alternative spellings; supply one.
Linear inequalities use `list(ineqA = A, ineqB = B)` for
$`A\theta + B \geq 0`$. Columns follow the full optimization parameter
order, including fixed coordinates. The stats constrained optimizer
requires an interior start and does not support equality constraints.

A `valid = function(param)` callback can reject invalid parameter
proposals by returning `FALSE`. Its errors propagate; it should be
deterministic.

## Residual parameters and identification

Residual parameters are fixed by default. Set `estimate_residual = TRUE`
to include them in estimation. This requires an identifiable model.

For example, with zero excitation and exponential residuals of rate
$`r`$, the event rate depends on the product $`r\mu`$. Estimating both
freely cannot separate them. Fixing the residual rate or an appropriate
kernel scale resolves that particular ambiguity. The package does not
automatically normalize an arbitrary residual law or diagnose every form
of nonidentification.

## Inspect the result

``` r

summary(fit)
#> FRPP maximum likelihood fit (stats, BFGS, bridge)
#> Log likelihood: 0.2958369
#> Convergence code: 0 (Normal convergence.)
#>        Estimate Std. Error z value Pr(>|z|)  
#> mu1      3.0000     1.7321  1.7321  0.08326 .
#> alpha1   0.0000     0.0000      NA       NA  
#> beta1    1.0000     0.0000      NA       NA  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
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
vcov(fit)
#>        mu1 alpha1 beta1
#> mu1      3      0     0
#> alpha1   0      0     0
#> beta1    0      0     0
c(AIC = AIC(fit), BIC = BIC(fit), events = nobs(fit))
#>       AIC       BIC    events 
#> 1.4083263 0.5069386 3.0000000
```

[`nobs()`](https://rdrr.io/r/stats/nobs.html) excludes the initial row.
AIC and BIC describe the fitted likelihood; they do not establish that
the model or its standard errors are appropriate.

[`fdiagnostics()`](https://ksublee.github.io/frpp/reference/ffit.md)
reports convergence and Wald covariance availability. Covariance is
withheld on optimizer failure, unavailable or nonpositive curvature, and
detected constraint boundaries. A reported `available_unvalidated`
covariance has not been checked for finite-sample confidence-interval
coverage.

Weak identification, long-tailed estimates and nonsmooth residual
densities can make Wald intervals unreliable. The moving transition in
`tzexp` deserves particular care. Profile likelihood and bootstrap
methods need separate validation; they are not implemented by
[`ffit()`](https://ksublee.github.io/frpp/reference/ffit.md).

## Restore and reuse a fitted model

``` r

estimated <- fitted_fspec(fit)
initial <- finitial_state(fit)
replayed <- do.call(logLik, c(list(object = estimated, inter_arrival = dt), initial))
c(optimized = as.numeric(logLik(fit)), replayed = as.numeric(replayed))
#> optimized  replayed 
#> 0.2958369 0.2958369
```

`fitted(fit)` also returns the estimated model.
[`finitial_state()`](https://ksublee.github.io/frpp/reference/ffit.md)
describes the start of the fitted observation window, not the state
after the final event. Continuing a path requires its terminal state.

Save results with [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) and
restore them with [`readRDS()`](https://rdrr.io/r/base/readRDS.html).
Captured callback environments remain references to external state;
saving a fit does not automatically freeze that external state.

## Validation of the migrated example

See the [10,000-row
comparison](https://ksublee.github.io/frpp/articles/legacy-fit-validation.md)
for matching `hfit()` and
[`ffit()`](https://ksublee.github.io/frpp/reference/ffit.md) estimates,
saved-fit replay and the unavailable Wald covariance in the
trapezoid-exponential example.
