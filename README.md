# frpp

An R package for flexible residual point processes, combining exponential
kernels from [emhawkes](https://github.com/ksublee/emhawkes) with custom residual distributions.

Use `fspec()` to define a model, `fsim()` to simulate, `logLik()` to evaluate
its likelihood, and `ffit()` to estimate parameters. This package is under development.

[Documentation](https://ksublee.github.io/frpp/)

## Installation

```r
install.packages("remotes")
remotes::install_github("ksublee/frpp", build_vignettes = FALSE)
```

## Example

```r
library(frpp)

kernel <- methods::new("hspec", mu = 0.3, alpha = 0.2, beta = 1)
law <- residual_law(c(shape = 2, rate = 2),
  stats::dgamma, stats::pgamma, stats::qgamma, stats::rgamma)
model <- fspec(kernel, law)
path <- fsim(model, size = 1000, seed = 73, lambda_component0 = 0)
```

The default `bridge` uses emhawkes 1.0.0 or later. Optional `state` simulation
and likelihood use the public kernel APIs in current GitHub development emhawkes.
See the documentation for support and timing limits.
