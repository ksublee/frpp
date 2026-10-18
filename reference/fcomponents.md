# Access Components of an FRPP Specification

Read or replace the kernel and residual law of an S4 `fspec`. Legacy S3
specifications are converted before access. Replacement validates the
whole specification without executing its callbacks.

## Usage

``` r
fkernel(object)
fkernel(object) <- value
fresidual(object)
fresidual(object) <- value
```

## Arguments

- object:

  An `fspec` or a compatible legacy FRPP specification.

- value:

  For `fkernel`, an `hspec` without residual callbacks. For `fresidual`,
  a law constructed by `residual_law`.

## Value

The requested component. Replacement returns an updated S4
specification.

## See also

[`fspec`](https://ksublee.github.io/frpp/reference/frpp_spec.md),
[`residual_law`](https://ksublee.github.io/frpp/reference/residual_law.md)

## Examples

``` r
m <- fspec(methods::new("hspec", mu = .4, alpha = .1, beta = 1),
  residual_law(c(rate = 1), stats::dexp, stats::pexp, stats::qexp, stats::rexp))
fkernel(m)
#> An object of class "hspec" of 1-dimensional Hawkes process
#> 
#> Slot mu: 
#>      [,1]
#> [1,]  0.4
#> 
#> Slot alpha: 
#>      [,1]
#> [1,]  0.1
#> 
#> Slot beta: 
#>      [,1]
#> [1,]    1
#> 
fresidual(m)$name
#> [1] "custom"
fresidual(m) <- residual_law(c(shape = 2, rate = 2), stats::dgamma,
  stats::pgamma, stats::qgamma, stats::rgamma)
methods::validObject(m)
#> [1] TRUE
```
